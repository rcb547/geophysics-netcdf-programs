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
#include <fstream>
using namespace netCDF;
using namespace netCDF::exceptions;

#define _PROGRAM_ "aseggdf2netcdf"
#define _VERSION_ "1.0"

#include "general_utils.h"
#include "file_utils.h"
#include "string_utils.h"
#include "blocklanguage.h"
#include "asciicolumnfile.h"
#include "streamredirecter.h"

#include "csvfile.h"
#include "geophysics_netcdf.hpp"

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

public:
	
	
	cASEGGDF2Converter(const std::string& datpath, const std::string& ncpath, const std::string& commandline){
		_GSTITEM_		
				
		DatPath = datpath;
		DatName = extractfilename_noextension(DatPath);
		DfnPath = extractfiledirectory(DatPath) + DatName + ".dfn";
		NCPath  = ncpath;
				
		std::string LogPath = NCPath + ".log";
		std::string WLogPath = NCPath + ".warn.log";
		
		std::ofstream wlog(WLogPath, std::ios::ate);
		cStreamRedirecter R(wlog,std::cerr);

		double t1 = gettime();
		glog.open(LogPath);	
		glog.logmsg("Program %s starting at %s\n", _PROGRAM_, timestamp().c_str());
		glog.logmsg("Version %s Compiled at %s on %s\n", _VERSION_, __TIME__, __DATE__);
		glog.logmsg("%s\n", commandline.c_str());
		glog.logmsg("Working directory: %s\n", getcurrentdirectory().c_str());
		convert_aseggdf2_file();		
		double t2 = gettime();
		glog.logmsg("Elapsed time = %.2lf\n", t2 - t1);
		glog.close();
	};

	~cASEGGDF2Converter(){						
		glog.logmsg("Finished at %s\n", timestamp().c_str());
		glog.close();					
	};

	bool convert_aseggdf2_file() {
		_GSTPUSH_

		if (exists(DatPath) == false) {
			glog.logmsg("Error 1: Data file %s does not exist\n", DatPath.c_str());
			return false;
		}

		if (exists(DatPath) == false) {
			glog.logmsg("Error 2: DFN file %s does not exist\n", DfnPath.c_str());
			return false;
		}

		glog.logmsg("Opening data file %s\n",DatPath.c_str());
		cAsciiColumnFile AF(DatPath);

		glog.logmsg("Parsing ASEGGDF2 header\n");
		AF.parse_dfn_header(DfnPath);

		//Force change attribute name "desc" or "DESC" to "description"
		//for (size_t i = 0; i < AF.fields.size(); i++) {
		//	for (size_t j = 0; j < AF.fields[i].atts.size(); j++) {
		//		if (tolower(AF.fields[i].atts[j].first) == "desc") {
		//			AF.fields[i].atts[j].first = "description";
		//		}
		//	}
		//}
						
		std::vector<unsigned int> line_number;
		std::vector<unsigned int> line_index_start;
		std::vector<unsigned int> line_index_count;

		int line_field_index = -1;
		std::string line_field_name = "";
		std::vector<std::string> cand = { "line", "linenumber", "line_number", "flightline", "fltline" };
		glog.logmsg("Determining line field name\n");
		for (size_t i = 0; i < cand.size(); i++) {
			glog.logmsg("\ttrying %s\n",cand[i].c_str());
			line_field_index = AF.fieldindexbyname(cand[i]);
			if (line_field_index >= 0) {
				line_field_name = AF.fields[line_field_index].name;
				break;
			}
		}
		
		if (line_field_index < 0) {
			std::string msg;
			msg += strprint("Could not determin the line number field'\n");						
			glog.errormsg(_SRC_+msg);
		}
		else {
			glog.logmsg("Using %s as the 'line number' field\n", line_field_name.c_str());
		}
		
		glog.logmsg("Scanning for line index\n");
		size_t npoints = AF.scan_for_line_index(line_field_index, line_index_start, line_index_count, line_number);
		glog.logmsg("Total number of points is %d\n", (int)npoints);
		glog.logmsg("Total number of lines is %d\n", (int)line_index_start.size());

		glog.logmsg("Scanning for groupby fields\n");
		std::vector<bool> isgroupby = AF.scan_for_groupby_fields(line_index_count);

		bool status = exists(extractfiledirectory(NCPath));
		if (status==false){
			makedirectorydeep(extractfiledirectory(NCPath));
		}
		
		glog.logmsg("Creating NetCDF file %s\n",NCPath.c_str());
		cGeophysicsNcFile ncFile(NCPath, NcFile::replace);

		glog.logmsg("Adding line index variables\n");
		ncFile.InitialiseNew(line_number, line_index_count);
				
		//Pre process the fields
		glog.logmsg("Pre processing fields\n");
		std::vector<nc_type> vartypes(AF.fields.size());
		std::vector<std::string> varnames(AF.fields.size());
		bool reported_nameswap = false;
		for (size_t fi = 0; fi < AF.fields.size(); fi++){
			cAsciiColumnField& f = AF.fields[fi];
										
			if (f.name == line_field_name) {
				std::string msg;
				msg += strprint("Skip processing field %3zu - %s (already in index)\n", fi, f.name.c_str());
				std::cout << msg << std::endl;
				glog.logmsg(msg);
				continue;
			}
			
			if (tolower(f.name) == "rt") {
				std::string msg;
				msg += strprint("Skip processing field %3zu - %s (is a RECORD_TYPR field)\n", fi, f.name.c_str());
				std::cout << msg << std::endl;
				glog.logmsg(msg);
				continue;
			}

			if (tolower(f.name) == "fltline") {
				std::string msg;
				msg += strprint("Skip processing field %3zu - %s (is the GEOSOFT fltline field)\n", fi, f.name.c_str());
				std::cout << msg << std::endl;
				glog.logmsg(msg);
				continue;
			}

			std::string longname = f.longname();						
			if (longname.size() == 0) {
				varnames[fi] = f.name;
			}
			else {
				varnames[fi] = longname;
				if (reported_nameswap == false) {
					std::string msg;
					msg += strprint("Warning: There are NAME=value pairs in the .dfn file'\n");
					msg += strprint("\tUsing those NAMES instead of the label immediately after the RT=;\n");
					std::cerr << msg << std::endl;
					reported_nameswap = true;					
				}
			}

			static const char space = ' ';
			size_t found = varnames[fi].find_first_of(space);
			if (found != std::string::npos) {
				std::string msg;
				msg += strprint("Error: There are spaces are in the variable name '%s'\n", varnames[fi].c_str());
				msg += strprint("\t The value 'NAME=value' pair in the .dfn file has got spaces in it\n");
				msg += strprint("\t If the NAMES are really descriptions, change the NAME=... to DESC=... instead\n");
				std::cerr << msg << std::endl;
				glog.errormsg(_SRC_+msg);
			}

			std::string& fieldname = varnames[fi];															
			size_t nbands = f.nbands;
			std::vector<NcDim> vardims;
			if (nbands > 1){								
				std::string dimname = f.get_att("second_dimension_name");
				if (dimname.size()==0) {
					std::string msg = strprint("Error 1: Multiband field %s has no SECOND_DIMENSION_NAME=... tag in .dfn file\n", varnames[fi].c_str());										
					glog.errormsg(_SRC_+msg);
				}
				NcDim dimband = ncFile.getDim(dimname);
				if (dimband.isNull()){
					std::vector<unsigned int> b = increment((unsigned int)nbands, (unsigned int)0, (unsigned int)1);
					dimband = ncFile.addDimVar(dimname, b);
				}
				dimband.isNull();
				vardims.push_back(dimband);
			}
			
			if (f.isinteger()){
				vartypes[fi] = NC_INT;
			}
			else if(f.isreal()){
				if (f.width > 8) {
					vartypes[fi] = NC_DOUBLE;
				}
				else {
					vartypes[fi] = NC_FLOAT;
				}				
			}
			else {
				std::string msg = strprint("Error unknown field datatype for %s\n", fieldname.c_str());				
				std::cerr << msg << std::endl;
				glog.errormsg(_SRC_+msg);
			}

			cGeophysicsVar gv;
			if (isgroupby[fi]) {
				gv = ncFile.addLineVar(varnames[fi], vartypes[fi], vardims);				
			}
			else {
				gv = ncFile.addSampleVar(varnames[fi], vartypes[fi], vardims);				
			}
			gv.add_original_dataset_fieldname(fieldname);
			
			for (const auto& [key, value] : f.atts) {				
				if (tolower(key) != "null") {
					//null is not relevant in the netcdf file
					gv.add_attribute(key, value);
				}
			}			
			
			std::string istr = "point";
			if(isgroupby[fi]) istr = "line";

			std::string tname = NcType(vartypes[fi]).getTypeClassName();			
			glog.logmsg("field index:%zu name:%s datatype:%s bands:%zu indexing:%s units:%s\n", fi + 1, fieldname.c_str(), tname.c_str(),nbands,istr.c_str(),f.units().c_str());
		}		

		size_t fi_line = AF.fieldindexbyname("line");				

		std::vector<std::vector<int>>    intfields;
		std::vector<std::vector<double>> dblfields;
		size_t lineindex = 0;
		AF.rewind();		
		AF.clear_currentrecord();
		glog.logmsg("Processing lines\n");
		size_t nsamples;
		while ((nsamples = AF.readnextgroup(fi_line, intfields, dblfields))){
			if (line_index_count[lineindex] != nsamples) {				
				std::string msg;
				msg += strprint("Error: number of samples read in from line does not match the index\n");
				msg += strprint("\tindex: %d and read in: %d\n", line_index_count[lineindex], nsamples);
				std::cerr << msg << std::endl;
				glog.errormsg(_SRC_+msg);
			}

			glog.logmsg("Processing line index:%zu linenumber:%u\n",lineindex+1,line_number[lineindex]);
			for (size_t fi = 0; fi < AF.fields.size(); fi++){				
				cAsciiColumnField& f = AF.fields[fi];
				std::string& vname = varnames[fi];				
				//std::cout << f.name << std::endl;								
				
				if (f.name == line_field_name) {
					continue;
				}

				if (f.ischar()) {
					continue;
				}

				if (tolower(f.name) == "rt") {					
					continue;
				}

				if (tolower(f.name) == "fltline") {					
					continue;
				}


				size_t nbands = AF.fields[fi].nbands;

				NcVar var = ncFile.getSampleVar(vname);
					
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
															
					if (AF.fields[fi].isinteger()){
						std::vector<int> data(nactive);
						int mv;
						cGeophysicsVar gv(&ncFile, var);
						mv = gv.missingvalue(mv);
						for (size_t si = 0; si < nactive; si++){							
							int& val = data[si] = intfields[fi][si*nbands + bi];
							if (!isdefined(val)) {
								val = mv;
							}
							else if (val == AF.fields[fi].nullvalue<int>()){
								val = mv;
							}
						}						
						var.putVar(startp, countp, data.data());						
					}
					else{						
						double mv;
						cGeophysicsVar gv(&ncFile,var);
						mv = gv.missingvalue(mv);												
						std::vector<double> data(nsamples);
						for (size_t si = 0; si < nactive; si++){							
							double& val = data[si] = dblfields[fi][si*nbands + bi];
							if (!isdefined(val)) {
								val = mv;
							}
							else if (val == AF.fields[fi].nullvalue<double>()){
								val = mv; 
							}
						}						
						var.putVar(startp, countp, data.data());						
					}
				}				
			}			
			lineindex++;
		}		
		add_global_attributes(ncFile);		
		glog.logmsg( "Conversion complete\n");
		_GSTPOP_
		return true;
	}
	
	bool add_global_attributes(cGeophysicsNcFile& ncFile) {
		_GSTITEM_
		ncFile.putAtt("CreationTime", timestamp());
		ncFile.putAtt("CreationMethod", "aseggdf2netcdf.exe");
		ncFile.putAtt("ASEGGDF2SourceDataFile", DatPath);
		ncFile.putAtt("ASEGGDF2SourceDFNFile", DfnPath);
		return true;
	}	
};

int main(int argc, char** argv)
{		
	_GSTITEM_		
	try
	{
		if (argc == 3) {			
			std::string datpath = argv[1];
			std::string ncpath  = argv[2];
			std::string cmdl = commandlinestring(argc, argv);
			cASEGGDF2Converter C(datpath,ncpath,cmdl);
			return 0;
		}
		else {
			std::cout << "Usage: " << extractfilename(argv[0]) << " datpath ncpath" << std::endl;
			return 1;
		}		
	}
	catch (NcException& e){
		_GSTPRINT_		
		std::cerr << e.what();
		//glog.logmsg(e.what());
		//glog.logmsg("\n");
		return 1;
	}
	catch (const std::runtime_error& e) {
		_GSTPRINT_
		std::cerr << e.what();
		//glog.logmsg(e.what());
		//glog.logmsg("\n");
		return 1;
	}
	catch (std::exception& e){
		_GSTPRINT_
		std::cerr << e.what();
		//glog.logmsg(e.what());
		//glog.logmsg("\n");
		return 1;
	}				
	return 0;
}





