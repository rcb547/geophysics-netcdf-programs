/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include <cstdio>
#include <netcdf>
#include <vector>
#include <limits>

using namespace netCDF;
using namespace netCDF::exceptions;

#ifdef USEGLOBALSTACKTRACE
#include "stacktrace.h"
class cStackTrace gtrace;
#endif

#define _PROGRAM_ "geophysicsnc2shape"
#define _VERSION_ "1.0"

#include "general_utils.h"
#include "file_utils.h"
#include "logger.h"
#include "ogr_utils.h"
#include "gdal_utils.h"
#include "geophysics_netcdf.hpp"

#ifdef HAVE_GDAL
	#include "crs.h"
#endif

class cLogger glog; //The instance of the global log file manager

class cNcToShapefileConverter {	
	std::string NCPath;
	std::string ShapePath;	
public:

	cNcToShapefileConverter(const std::string& ncfilepath, const std::string& shapefilepath) {
		_GSTITEM_;
		NCPath    = fixseparator(ncfilepath);
		ShapePath = fixseparator(shapefilepath);						
		bool status = process();	
		if (status == false) {
			glog.logmsg("Error 0: creating shapefile %s from %s\n",ShapePath.c_str(),NCPath.c_str());
		}
		glog.close();			
	};

	~cNcToShapefileConverter() {};

	std::string look_for_var(cGeophysicsNcFile& N,  const std::vector<std::string>& candidates)
	{
		std::vector<NcVar> v = N.getAllVars();		
		for(size_t ci=0; ci<candidates.size(); ci++){
			for (size_t vi = 0; vi < v.size(); vi++) {
				if(strcasecmp(v[vi].getName(), candidates[ci]) == 0){
					return v[vi].getName();
				}
			}
		}
		return std::string();
	}

	std::vector<double> get_linetype(cGeophysicsNcFile& N)
	{		
		std::vector<std::string> candidates;
		candidates.push_back("linetype");
		candidates.push_back("line_type");
		candidates.push_back("ltype");
		std::string vname = look_for_var(N, candidates);
		std::vector<double> ltype;

		if (vname.size() > 0) {
			NcVar var = N.getVar(vname);
			if (N.isLineVar(var)) {
				cLineVar v = N.getLineVar(vname);
				v.getAll(ltype);
			}
			else if (N.isSampleVar(var)) {
				ltype.resize(N.nlines());
				cSampleVar v = N.getSampleVar(vname);
				for (size_t li = 0; li < N.nlines(); li++) {
					v.getSample(li, 0, 0, ltype[li]);
				}
			}
		}
		return ltype;
	}

	bool process() {	
		if (!exists(extractfiledirectory(ShapePath))){
			makedirectorydeep(extractfiledirectory(ShapePath));
		}
		cGeoDataset D = cGeoDataset::create_shapefile(ShapePath);
		cLayer L = D.create_layer("flight_lines", OGRwkbGeometryType::wkbLineString);
		std::vector<cAttribute> atts;
		atts.push_back(cAttribute("linenumber", (int)0));
		atts.push_back(cAttribute("linetype", (int)0));
		L.add_fields(atts);

		cGeophysicsNcFile N(NCPath);
		std::vector<unsigned int> ln;
		N.getLineNumbers(ln);
		const size_t nl = N.nlines();
				
		std::string xvarname;				
		std::vector<std::string> xcand = { "longitude","longitude_gda94" };
		for (size_t i = 0; i < xcand.size(); i++) {
			xvarname = xcand[i];
			if (N.hasVarCaseInsensitive(xcand[i])) {				
				xvarname = xcand[i];
				break;
			}
		}
		if(xvarname.size()==0){				
			std::string msg = _SRC_ + strprint("Could not find field longitude or longitude_gda94 (%s)\n", NCPath.c_str());
			throw(std::runtime_error(msg));
		}

		std::string yvarname;
		std::vector<std::string> ycand = { "latitude","latitude_gda94" };
		for (size_t i = 0; i < ycand.size(); i++) {
			yvarname = ycand[i];
			if (N.hasVarCaseInsensitive(ycand[i])){				
				yvarname = ycand[i];
				break;
			}
		}
		if (yvarname.size() == 0) {
			std::string msg = _SRC_ + strprint("Could not find field latitude or latitude_gda94 (%s)\n", NCPath.c_str());
			throw(std::runtime_error(msg));
		}

		std::vector<double> ltype = get_linetype(N);				
		std::vector<double> lkm0(nl);
		std::vector<double> lkm2(nl);
		std::vector<double> lkm4(nl);
		for (size_t li = 0; li < nl; li++) {
			lkm0[li] = 0.0;
			lkm2[li] = 0.0;
			lkm4[li] = 0.0;

			std::vector<double> x;
			std::vector<double> y;
			
			N.getDataByLineIndex(xvarname, li, x);			
			N.getDataByLineIndex(yvarname, li, y);

			std::vector<double> xout;
			std::vector<double> yout;

			double null = defaultmissingvalue(ncDouble);

			int ns = (int) x.size();
			int k = 0;
			int kstart, kend;
			while (x[k] == null || y[k] == null) {				
				k++;
				if (k == ns)break;
			}
			kstart = k;
			
			k = ns - 1;
			while (x[k] == null || y[k] == null) {				
				k--;
				if (k == -1)break;
			}
			kend = k;
			
			if(kend >= kstart){
				int minpoints = 20;
				int ss = (kend - kstart) / (minpoints - 2);
				if (ss < 1) ss = 1;
				for (k = kstart; k < kend; k += ss) {
					if(x[k] != null && y[k] != null) {
						xout.push_back(x[k]);
						yout.push_back(y[k]);						
					}					
				}
				xout.push_back(x[kend]);
				yout.push_back(y[kend]);			

				atts[0].value = (int)ln[li];
				atts[1].value = (int)0;
				if (ltype.size() == nl) {
					atts[1].value = (int)ltype[li];
				}
				L.add_linestring_feature(atts, xout, yout);
			}
		}				
		return true;
	}
};

int main(int argc, char** argv)
{
	_GSTITEM_

	GDALAllRegister();

	glog.logmsg(0,"Program %s \n", _PROGRAM_);
	glog.logmsg(0,"Version %s Compiled at %s on %s\n", _VERSION_, __TIME__, __DATE__);
	glog.logmsg(0,"Working directory %s\n", getcurrentdirectory().c_str());

	try
	{		
		if (argc == 3) {
			std::string NCPath    = argv[1];
			std::string ShapePath = argv[2];									
			std::cout << NCPath << " " << ShapePath << std::endl << std::flush;
			cNcToShapefileConverter C(NCPath, ShapePath);			
			glog.logmsg(0, "Finished\n");
		}
		else if (argc == 4) {			
			std::string ncdir = argv[1];
			std::string shapedir = argv[2];
			std::string listfile = argv[3];
			std::ifstream file(listfile);
			addtrailingseparator(ncdir);
			addtrailingseparator(shapedir);
			int k = 0;
			while (file.eof() == false) {
				std::string s;
				file >> s;
				s = trim(s);
				if (s.size() > 0 && s[0] != '#') {
					sFilePathParts fpp = getfilepathparts(s);
					std::string NCPath = ncdir + fpp.directory + fpp.prefix + ".nc";
					std::string ShapePath = shapedir + fpp.directory + fpp.prefix + ".shp";
					std::cout << NCPath << " " << ShapePath << std::endl << std::flush;
					cNcToShapefileConverter C(NCPath, ShapePath);
					k++;
				}
			}
			glog.logmsg(0, "Finished\n");
		}
		else{
			std::cout << "Usage: " << extractfilename(argv[0]) << " ncfile shapefile" << std::endl;
			std::cout << "   or: " << extractfilename(argv[0]) << " ncfiles_directory shapefiles_directory list_of_ncfiles.txt" << std::endl;
		}
	}
	catch (NcException& e)
	{
		_GSTPRINT_
		std::cout << e.what() << std::endl;		
		return 1;
	}
	catch (std::exception& e)
	{
		_GSTPRINT_
		std::cout << e.what() << std::endl;		
		return 1;
	}		
	return 0;
}



