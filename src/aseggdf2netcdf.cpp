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
#include <algorithm>
using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;

#include <fstream>
#include "general_utils.h"
#include "file_utils.h"
#include "asciicolumnfile.h"
#include "gdal_utils.h"

#include "metadata.h"
#include "crs.h"
#include "csvfile.h"
#include "geophysics_netcdf.h"

#include "stacktrace.h"
#ifdef USEGLOBALSTACKTRACE	
cStackTrace gtrace;
#endif

#ifdef USEGLOBALSTACKTRACE
	#include "stacktrace.h"
	cStackTrace globalstacktrace;
#endif
class cLogger glog; //The instance of the global log file manager

class cASEGGDF2Converter{
	
	std::string DatPath;
	std::string DatName;
	std::string DfnPath;	
	std::string NCPath;
	std::string LogFile="aseggdf2netcdf.log";

public:
	
	
	cASEGGDF2Converter(const std::string& datpath, const std::string& ncpath, const std::string& controlfile){
		_GSTITEM_		
		cBlock B(controlfile);
		
		DatPath = datpath;
		DatName = extractfilename_noextension(DatPath);
		DfnPath = extractfiledirectory(DatPath) + DatName + ".dfn";
		NCPath  = ncpath;

		std::string logfile = B.getstringvalue("LogFile");
		
		if (isdefined(logfile)) {
			LogFile = logfile;
		}
		glog.open(LogFile);				
		glog.logmsg("Program starting at %s\n", timestamp().c_str());
		glog.log(B.get_as_string());

		convert_aseggdf2_file();				
		glog.close();
	};

	~cASEGGDF2Converter(){						
		glog.logmsg("Finished at %s\n", timestamp().c_str());
		glog.close();					
	};

	size_t scan_for_line_index(const cAsciiColumnFile& D, std::vector<unsigned int>& line_index_start, std::vector<unsigned int>& line_index_count, std::vector<unsigned int>& line_number)
	{
		_GSTITEM_
		size_t fi = D.fieldindexbyname("line");
		size_t i1 = D.fields[fi].startchar;
		size_t i2 = D.fields[fi].endchar;

		unsigned int n = 0;
		std::ifstream in(DatPath);
		std::string s;

		int width = (int)(i2 - i1 + 1);
		unsigned int lastline = 0;
		size_t nlines = 0;
		while (std::getline(in, s)){
			std::string t = s.substr(i1, width);
			unsigned int lnum = atoi(t.data());
			if (lnum != lastline){
				line_number.push_back(lnum);
				line_index_start.push_back(n);
				line_index_count.push_back(1);
				lastline = lnum;
			}
			else{
				line_index_count.back()++;
			}
			n++;
		}
		return n;
	}

	std::vector<bool>  scan_for_groupby_fields(const cAsciiColumnFile& D, const std::vector<unsigned int>& line_index_count)
	{		
		_GSTITEM_
		std::vector<bool> groupby(D.fields.size(),true);
		std::ifstream infile(DatPath);
		std::string s, t;			
		int nl = std::min((int)2,(int)line_index_count.size());
		for (size_t li = 0; li < nl; li++){						
			std::getline(infile, s);
			for (size_t si = 1; si < line_index_count[li]; si++){				
				std::getline(infile, t);
				for (size_t fi = 0; fi < D.fields.size(); fi++){					
					if (groupby[fi] == true){						
						size_t i1 = D.fields[fi].startchar;
						size_t i2 = D.fields[fi].endchar;
						size_t width = i2 - i1 + 1;
						std::string a = s.substr(i1, width);
						std::string b = t.substr(i1, width);
						if (a != b) groupby[fi] = false;						
					}
				}
			}
		}			
		return groupby;		
	}
	
	bool convert_aseggdf2_file() {
		_GSTITEM_

			if (exists(DatPath) == false) {
				glog.logmsg("Error 1: Data file %s does not exist\n", DatPath.c_str());
				return false;
			}

		if (exists(DatPath) == false) {
			glog.logmsg("Error 2: DFN file %s does not exist\n", DfnPath.c_str());
			return false;
		}

		glog.logmsg("Opening data file %s\n",DatPath.c_str());
		cAsciiColumnFile D(DatPath);

		glog.logmsg("Parsing ASEGGDF2 header\n");
		D.parse_aseggdf2_header(DfnPath);

		std::vector<unsigned int> line_number;
		std::vector<unsigned int> line_index_start;
		std::vector<unsigned int> line_index_count;
		glog.logmsg("Scanning for line index\n");
		size_t npoints = scan_for_line_index(D, line_index_start, line_index_count, line_number);
		glog.logmsg("Total number of points is %lu\n", npoints);

		glog.logmsg("Scanning for groupby fields\n");
		std::vector<bool> isgroupby = scan_for_groupby_fields(D, line_index_count);

		bool status = exists(extractfiledirectory(NCPath));
		if (status==false){
			makedirectorydeep(extractfiledirectory(NCPath));
		}
		
		glog.logmsg("Creating NetCDF file %s\n",NCPath.c_str());
		cGeophysicsNcFile ncFile(NCPath, NcFile::replace);

		glog.logmsg("Adding line index variables\n");
		ncFile.InitialiseNew(line_number, line_index_count);
						
		//Pre process the fields
		for (size_t fi = 0; fi < D.fields.size(); fi++){
			cAsciiColumnField& f = D.fields[fi];
			std::string fieldname = f.name;
			if (fieldname == DN_LINE || fieldname == "Line"){
				glog.logmsg("Skip processing field %3zu - %s (already in index)\n", fi, fieldname.c_str());
				continue;
			}			
			else{
				glog.logmsg("Pre processing field  %3lu - %s\n", fi, fieldname.c_str());
			}
						
			size_t nbands = f.nbands;

			std::vector<NcDim> vardims;
			if (nbands > 1){
				std::string dimname = "window";
				NcDim dimband = ncFile.getDim(dimname);
				if (dimband.isNull()){
					std::vector<unsigned int> b = increment((unsigned int)nbands, (unsigned int)0, (unsigned int)1);
					dimband = ncFile.addDimVar(dimname, b);
				}
				bool status = dimband.isNull();
				vardims.push_back(dimband);
			}

			std::string varname = fieldname;
			nc_type vartype;

			if (f.isinteger()){
				vartype = NC_INT;
			}
			else if(f.isreal()){
				if (f.fmtwidth > 8) {
				    vartype = NC_DOUBLE;
				}
				else {
					vartype = NC_FLOAT;
				}				
			}
			else {
				glog.logmsg("Error unknown field datatype for %s\n", fieldname.c_str());
			}

			if (isgroupby[fi]){
				cLineVar var = ncFile.addLineVar(varname, vartype, vardims);				
				//var.add_long_name(fieldname);
				var.add_original_name(fieldname);
				var.add_units(f.units);
			}
			else{
				cSampleVar var = ncFile.addSampleVar(varname, vartype, vardims);
				//var.add_long_name(fieldname);
				var.add_original_name(fieldname);
				var.add_units(f.units);

				//if (strcasecmp(fieldname.c_str(), "latitude") == 0) var.add_standard_name("latitude");
				//else if (strcasecmp(fieldname.c_str(), "Latitude") == 0) var.add_standard_name("latitude");
				//else if (strcasecmp(fieldname.c_str(), "longitude") == 0) var.add_standard_name("longitude");
				//else if (strcasecmp(fieldname.c_str(), "Longitude") == 0) var.add_standard_name("longitude");
				//else;				
			}
						
			NcType t(vartype);
			std::string tname = NcType(vartype).getTypeClassName();
			glog.logmsg("type %s\n", tname.c_str());
			
		}

		size_t fi_line = D.fieldindexbyname("line");				

		std::vector<std::vector<int>>    intfields;
		std::vector<std::vector<double>> dblfields;
		size_t lineindex = 0;
		glog.logmsg("Processing lines\n");
		while (size_t nsamples = D.readnextgroup(fi_line, intfields, dblfields)){			
			for (size_t fi = 0; fi < D.fields.size(); fi++){				
				std::string fieldname = D.fields[fi].name;

				//if (fieldname == DN_LINE) continue;
				//if (fieldname == "Line") continue;
				if (fieldname == DN_LINE || fieldname == "Line") continue;

				size_t nbands = D.fields[fi].nbands;

				NcVar var = ncFile.getSampleVar(fieldname);
				
				std::vector<size_t> startp(2);
				std::vector<size_t> countp(2);								
				for (size_t bi = 0; bi < nbands; bi++){
					size_t nactive;
					if (isgroupby[fi]){
						nactive = 1;
						startp[0] = lineindex;
						startp[1] = bi;
						countp[0] = 1;
						countp[1] = 1;
					}
					else{
						nactive = nsamples;
						startp[0] = line_index_start[lineindex];
						startp[1] = bi;
						countp[0] = line_index_count[lineindex];
						countp[1] = 1;
					}
															
					if (D.fields[fi].isinteger()){
						std::vector<int> data(nactive);
						for (size_t si = 0; si < nactive; si++){							
							data[si] = intfields[fi][si*nbands + bi];
							if (data[si] == D.fields[fi].nullvalue){
								data[si] = defaultmissingvalue(ncInt);
							}
						}
						var.putVar(startp, countp, data.data());
					}
					else{
						std::vector<double> data(nsamples);
						for (size_t si = 0; si < nactive; si++){
							data[si] = dblfields[fi][si*nbands + bi];
							if (data[si] == D.fields[fi].nullvalue){
								data[si] = defaultmissingvalue(ncDouble); 
							}
						}
						var.putVar(startp, countp, data.data());
					}
				}				
			}			
			lineindex++;
		}		
		
		glog.logmsg( "Conversion complete\n");

		return true;
	}
	
	/*
	bool add_global_metadata(cGeophysicsNcFile& ncFile, const cMetaDataRecord& m){		
		_GSTITEM_
		for (size_t j = 0; j < m.header.size(); j++){
			if (m.header[j][0] == '/')continue;

			if (strcasecmp(m.header[j], "geoscience_australia_airborne_survey_project_number") == 0){
				int pnum = atoi(m.values[j].c_str());
				ncFile.putAtt(m.header[j].c_str(), ncInt, pnum);
			}
			else if (strcasecmp(m.header[j], "nominal_minimum_line_spacing") == 0){
				std::string s = m.values[j] + " m";
				ncFile.putAtt(m.header[j].c_str(), s);
			}
			else if (strcasecmp(m.header[j], "nominal_maximum_line_spacing") == 0){
				std::string s = m.values[j] + " m";
				ncFile.putAtt(m.header[j].c_str(), s);
			}
			else if (strcasecmp(m.header[j], "nominal_height_above_ground") == 0){
				std::string s = m.values[j] + " m";
				ncFile.putAtt(m.header[j].c_str(), s);
			}
			else if (strcasecmp(m.header[j], "nominal_height_above_sealevel") == 0){
				std::string s = m.values[j] + " m";
				ncFile.putAtt(m.header[j].c_str(), s);
			}
			else if (strcasecmp(m.header[j], "acquisition_start_date") == 0){
				std::string s = standardize_date(m.values[j]);
				ncFile.putAtt(m.header[j].c_str(), s);
			}
			else if (strcasecmp(m.header[j], "acquisition_end_date") == 0){
				std::string s = standardize_date(m.values[j]);
				ncFile.putAtt(m.header[j].c_str(), s);
			}
			else{
				ncFile.putAtt(m.header[j].c_str(), m.values[j].c_str());
			}
		}
		return true;
	}
	*/
};

int main(int argc, char** argv)
{		
	_GSTITEM_		
	try
	{
		if (argc == 4) {
			double t1 = gettime();			
			std::string datpath = argv[1];
			std::string ncpath  = argv[2];
			std::string controlfile = argv[3];
			cASEGGDF2Converter C(datpath,ncpath,controlfile);
			double t2 = gettime();
			printf("Reading data elapsed time = %.2lf\n", t2 - t1);
		}
		else {
			std::cout << "Usage: " << extractfilename(argv[0]) << " datpath ncpath control_file_name" << std::endl;
			return 1;
		}
	}
	catch (NcException& e)
	{
		_GSTPRINT_		
		glog.logmsg(e.what());
		glog.logmsg("\n");
		return 1;
	}
	catch (std::exception& e)
	{
		_GSTPRINT_
		glog.logmsg(e.what());
		glog.logmsg("\n");
		return 1;
	}			
	return 0;
}





