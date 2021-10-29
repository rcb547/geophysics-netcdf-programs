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

#define _PROGRAM_ "intrepid2netcdf"
#define _VERSION_ "1.0"

#include "general_utils.h"
#include "file_utils.h"
#include "blocklanguage.h"
#include "intrepid.h"
#include "streamredirecter.h"

#include "metadata.h"
#include "csvfile.h"
#include "logger.h"
#include "geophysics_netcdf.h"
#ifdef HAVE_GDAL
	#include "crs.h"
#endif

class cLogger glog; //The instance of the global log file manager

class cIntrepidToNetCDFConverter {	
	std::string IntrepiDatabasePath;
	std::string NCPath;
	bool OverWriteExistingNcFiles = true;
	
public:

	cIntrepidToNetCDFConverter(const std::string& intrepiddatabasepath, const std::string& ncfilepath, std::string& commandline) {
		_GSTITEM_;
		IntrepiDatabasePath = fixseparator(intrepiddatabasepath);
		NCPath = fixseparator(ncfilepath);
				
		std::string LogPath = NCPath + ".log";
		std::string WLogPath = NCPath + ".warn.log";

		std::ofstream wlog(WLogPath, std::ios::ate);
		cStreamRedirecter R(wlog, std::cerr);
		
		double t1 = gettime();
		glog.open(LogPath);
		glog.logmsg("Program %s starting at %s\n", _PROGRAM_, timestamp().c_str());		
		glog.logmsg("Version %s Compiled at %s on %s\n", _VERSION_, __TIME__, __DATE__);
		glog.logmsg("%s\n", commandline.c_str());
		glog.logmsg("Working directory: %s\n", getcurrentdirectory().c_str());		
		bool status = process();	
		if (status == false) {
			std::string msg = strprint("Error 0: converting %s to %s\n", intrepiddatabasepath.c_str(), ncfilepath.c_str()); 
			glog.logmsg(msg);	
			std::cerr << msg << std::endl;
		}
		double t2 = gettime();
		glog.logmsg("Elapsed time = %.2lf\n", t2 - t1);		
		glog.close();		
	};

	~cIntrepidToNetCDFConverter() {


	};

	bool process() {

		std::string IDBPath = ILDataset::dbdirpath(IntrepiDatabasePath);
		std::string IDBName = ILDataset::dbname(IntrepiDatabasePath);

		if (exists(NCPath)){
			if (OverWriteExistingNcFiles == false) {
				glog.logmsg("Warning 1: NetCDF file %s already exists - skipping this database\n", NCPath.c_str());
				return false;
			}
		}

		glog.logmsg("\nConverting database: %s\n", IntrepiDatabasePath.c_str());
		bool datasetexists = exists(IntrepiDatabasePath);
		if (datasetexists == false) {
			glog.logmsg("Error 1: Database does not exist: %s\n", IntrepiDatabasePath.c_str());
			return true;
		}

		glog.logmsg("\nOpening Intrepid database\n");
		ILDataset D(IntrepiDatabasePath);
		if (D.ispointdataset()) {
			glog.logmsg("Error 2: point databases not supported - skipping %s\n",IntrepiDatabasePath.c_str());
			return true;
		}
		else if(D.valid == false){
			glog.logmsg("Error 3: problem opening database - skipping %s\n",IntrepiDatabasePath.c_str());
			return true;
		}

		glog.logmsg("\nGetting the line numbers\n");
		std::string linenumberfieldname;
		D.getlinenumberfieldname(linenumberfieldname);
		if (D.fieldexists(linenumberfieldname) == false) {
			glog.logmsg("Error 4: could not determine the LineNumber field in the SurveyInfo file for - skipping %s\n", IntrepiDatabasePath.c_str());
			return true;
		}
		
		std::vector<size_t> linenumbers;
		if (D.getlinenumbers(linenumbers) == false) {					
			glog.logmsg("Error 5: could not determine the line numbers - skipping % s\n", IntrepiDatabasePath.c_str());
			return true;
		}					
		
		if (D.hassurveyinfokey_and_fieldexists("X") == false) {
			glog.logmsg("Warning 2: could not determine the X field in the SurveyInfo file\n");
		}

		if (D.hassurveyinfokey_and_fieldexists("Y") == false) {
			glog.logmsg("Warning 3: could not determine the Y field in the SurveyInfo file\n");
		}
		
		glog.logmsg("Creating NetCDF file: %s\n",NCPath.c_str());
		cGeophysicsNcFile ncFile(NCPath, NcFile::replace);

		glog.logmsg("\nAdding the line index variable\n");				
		std::vector<size_t> count = D.linesamplecount();									
		ncFile.InitialiseNew(linenumbers, count);			
		
		glog.logmsg("\nAdding global attributes\n");
		add_global_attributes(ncFile);

		glog.logmsg("\nAdding groupby varaibles\n");
		add_groupbyline_variables(ncFile, D);

		glog.logmsg("\nAdding indexed varaibles\n");
		add_indexed_variables(ncFile, D);

		glog.logmsg("\nConversion complete\n");
		return true;
	}

	NcType nc_datatype(const ILField& F)
	{
		_GSTITEM_
		if (F.getType().isubyte()) return NcType(ncUbyte);
		else if (F.getType().isshort()) return NcType(ncShort);
		else if (F.getType().isint()) return NcType(ncInt);
		else if (F.getType().isfloat())return NcType(ncFloat);
		else if (F.getType().isdouble())return NcType(ncDouble);
		else if (F.getType().isstring()) {
			return NcType(ncString);
		}
		else {
			std::string msg = strprint("Error 6: Unknown Intrepid data type in %s\n",F.datafilepath().c_str());
			glog.logmsg(msg);
			throw(msg);
		}
	}

	bool set_intrepid_nullvalue(cGeophysicsVar& v)
	{
		_GSTITEM_
		NcType t = v.getType();
		if (t == ncUbyte) v.add_missing_value(IDataType::ubytenull());
		else if (t == ncShort) v.add_missing_value(IDataType::shortnull());
		else if (t == ncInt) v.add_missing_value(IDataType::intnull());		
		else if (t == ncFloat) v.add_missing_value(IDataType::floatnull());
		else if (t == ncDouble) v.add_missing_value(IDataType::doublenull());
		else {
			std::string msg = strprint("Error 7: Unsupported data type %s for field %s in %s\n",t.getName().c_str(),v.getName().c_str(),IntrepiDatabasePath.c_str());
			glog.logmsg(msg);
			throw(msg);
		}
		return true;
	}

	bool add_groupbyline_variables(cGeophysicsNcFile& ncFile, ILDataset& D)
	{	
		_GSTITEM_
		if (D.valid == false)return false;
		size_t nlines = D.nlines();
		NcDim  dim_line = ncFile.getDim(DN_LINE);

		std::string linenumberfield;
		D.getlinenumberfieldname(linenumberfield);
		size_t fi = 0;
		for (auto it = D.Fields.begin(); it != D.Fields.end(); ++it) {
			fi++;
			ILField& F = *it;
			if (F.isgroupbyline() == false) continue;			
			if (F.getName() == linenumberfield) continue;
			if (F.getTypeId() == IDataType::ID::UNKNOWN) {
				glog.logmsg("Warning 4: skipping field %s: unsupported datatype\n", F.datafilepath().c_str());
				continue;
			}			
			
			glog.logmsg("Converting field %s\n", F.getName().c_str());
			std::vector<NcDim> dims;
			if (F.nbands() > 1) {
				std::string dimname = "nbands_" + F.getName();
				NcDim dim_band = ncFile.addDim(dimname, F.nbands());
				dims.push_back(dim_band);
			}

			std::vector<int> vstringasint;	
			nc_type outdatatype = nc_datatype(F).getId();
			if(F.getTypeId() == IDataType::ID::STRING){
				D.getgroupbydata(F, vstringasint);
				outdatatype = ncInt.getId();
				glog.logmsg("Warning 5: Converting field %s: with STRING datatype to 'int' datatype\n", F.datafilepath().c_str());
			}
			
			cLineVar var = ncFile.addLineVar(F.getName(), NcType(outdatatype), dims);
			for (size_t li = 0; li < nlines; li++) {
				ILSegment S(F,li);

				if (S.readbuffer() == false) {
					glog.logmsg("Error 8: could not read buffer for line sequence number %zu in field %s\n", li, F.getName().c_str());
					return false;
				}				
				change_fillvalues(S);

				std::vector<size_t> startp(2);
				std::vector<size_t> countp(2);

				//Point dimension
				startp[0] = li;
				countp[0] = 1;

				//Band dimension
				startp[1] = 0;
				countp[1] = S.nbands();
				if (F.getTypeId() == IDataType::ID::STRING) {
					var.putVar(startp, countp, (void*)&(vstringasint[li]));
				}
				else {
					var.putVar(startp, countp, S.pvoid_groupby());
				}
			}	
			add_field_attributes(F, var);
		}
		return true;
	}

	bool add_indexed_variables(cGeophysicsNcFile& ncFile, ILDataset& D)
	{
		_GSTITEM_
		if (D.valid == false)return false;
		size_t nlines = D.nlines();

		size_t fi = 0;
		for (auto it = D.Fields.begin(); it != D.Fields.end(); ++it) {
			fi++;
			ILField& F = *it;
			if (F.isgroupbyline() == true) continue;			
			
			if (ncFile.getVar(F.getName()).isNull() == false) {
				glog.logmsg("Warning 6: variable name %s already exists in this NC file - skipping field %s\n", F.getName().c_str(), F.datafilepath().c_str());
				return false;
			}

			if (F.getTypeId() == IDataType::ID::UNKNOWN) {
				glog.logmsg("Warning 4: skipping field %s: unsupported datatype\n", F.datafilepath().c_str());
				continue;
			}

			nc_type outdatatype = nc_datatype(F).getId();
			if (F.getTypeId() == IDataType::ID::STRING) {				
				outdatatype = ncInt.getId();
				glog.logmsg("Warning 5: Converting field %s: with STRING datatype to 'int' datatype\n", F.datafilepath().c_str());
			}
						
			glog.logmsg("Converting field %s\n", F.getName().c_str());
			std::vector<NcDim> dims;
			if (F.nbands() > 1) {
				std::string dimname = "nbands_" + F.getName();
				NcDim dim_band = ncFile.addDim(dimname, F.nbands());
				dims.push_back(dim_band);
			}
					   			 		  		  
			cSampleVar var = ncFile.addSampleVar(F.getName(), NcType(outdatatype), dims);
			size_t startindex = 0;
			for (size_t li = 0; li < nlines; li++) {
				ILSegment S(F,li);

				if (S.readbuffer() == false) {
					glog.logmsg("Error 9: could not read buffer for line sequence number %zu in field %s\n", li, F.datasetpath().c_str());
					return false;
				}
				change_fillvalues(S);
				
				std::vector<size_t> startp(2);
				std::vector<size_t> countp(2);

				//Sample dimension
				startp[0] = startindex;
				countp[0] = S.nsamples();

				//Band dimension
				startp[1] = 0;
				countp[1] = S.nbands();

				if (F.getTypeId() == IDataType::ID::STRING) {
					std::vector<int> vstringasint;
					S.getband(vstringasint, 0);
					var.putVar(startp, countp, (void*)vstringasint.data());
					
					//std::vector<std::string> svec;
					//bool stat = S.getband(svec, 0);
					//std::vector<const char*> p(S.nsamples());
					//for (size_t si = 0; si < S.nsamples(); si++) {
					//	p[si] = svec[si].c_str();
					//}
					//var.putVar(startp, countp, p.data());
				}
				else {
					var.putVar(startp, countp, S.pvoid());
				}
				startindex += S.nsamples();
			}
			add_field_attributes(F, var);
		}
		return true;
	}
	
	void change_fillvalues(ILSegment& S) {
		//Replace the nulls with NetCDF default fill values
		if (S.getType().isubyte()) {
			S.change_nullvalue(defaultmissingvalue(ncUbyte));
		}
		else if (S.getType().isshort()) {
			S.change_nullvalue(defaultmissingvalue(ncShort));
		}
		else if (S.getType().isint()) {
			S.change_nullvalue(defaultmissingvalue(ncInt));
		}
		else if (S.getType().isfloat()) {
			S.change_nullvalue(defaultmissingvalue(ncFloat));
		}
		else if (S.getType().isdouble()) {
			S.change_nullvalue(defaultmissingvalue(ncDouble));
		}
		else if (S.getType().isstring()) {
			S.change_nullvalue(defaultmissingvalue(ncString));
		}
		else {
			std::string msg = strprint("Error 10: Unsupported data type %s for field %s in %s\n",S.getType().getName().c_str(),S.getField().getName().c_str(),IntrepiDatabasePath.c_str());
			glog.logmsg(msg);
			throw(msg);
		}
	}

	void add_field_attributes(const ILField& F, cGeophysicsVar& var) {
		if (F.Datum.size() > 0) var.add_attribute("IntrepidDatumString", F.Datum);
		if (F.Projection.size() > 0) var.add_attribute("IntrepidProjectionString", F.Projection);
		if (F.CoordinateType.size() > 0) var.add_attribute("IntrepidCoordinateTypeString", F.CoordinateType);
		
		std::string aliaskey;
		if (F.getDataset().fieldalias(F.getName(), aliaskey)) {
			var.add_attribute("IntrepidAlias", aliaskey);
		}
	}

	bool add_global_attributes(cGeophysicsNcFile& ncFile) {
		_GSTITEM_
		ncFile.putAtt("CreationTime", timestamp());
		ncFile.putAtt("CreationMethod", "intrepid2netcdf.exe");		
		ncFile.putAtt("IntrepidSourceDataset", IntrepiDatabasePath);
		return true;
	}	
};

int main(int argc, char** argv)
{
	_GSTITEM_				
	std::string cmdl = commandlinestring(argc, argv);
	try
	{
		if (argc == 3) {
			std::string dbname = argv[1];
			std::string ncname = argv[2];
			cIntrepidToNetCDFConverter C(dbname, ncname, cmdl);
			glog.logmsg("Finished\n");
		}
		else if (argc == 4) {		
			std::string dbdir = argv[1];
			std::string ncdir = argv[2];
			std::string listfile = argv[3];
			std::ifstream file(listfile);
			addtrailingseparator(dbdir);
			addtrailingseparator(ncdir);
			while (file.eof() == false) {
				std::string db;
				file >> db;
				db = trim(db);
				if (db.size() > 0 && db[0] != '#') {
					sFilePathParts fpp = getfilepathparts(db);
					std::string dbname = dbdir + fpp.directory + fpp.prefix;
					std::string ncname = ncdir + fpp.directory + fpp.prefix + ".nc";					
					std::cout << dbname << " " << ncname << std::endl << std::flush;
					cIntrepidToNetCDFConverter C(dbname, ncname, cmdl);										
				}
			}			
		}
		else {			
			std::cout << "Usage: " << extractfilename(argv[0]) << " input_database output_ncfile" << std::endl;			
			std::cout << "   or: " << extractfilename(argv[0]) << " databases_dir ncfiles_dir list_of_databases.txt" << std::endl;
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


