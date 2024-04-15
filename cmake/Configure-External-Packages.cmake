message(STATUS "")
message(STATUS "== Configuring external packages =====================================")
list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Find PkgConfig
find_package(PkgConfig QUIET)

# Configure NetCDF C++ Libraries
if(${WITH_NETCDF})
	message(STATUS "\nChecking for NETCDFCXX")
	include(cmake/Configure-NETCDF_CXX.cmake)
	if(NETCDFCXX_FOUND)
		message(STATUS "NETCDFCXX was found")
	else()
		message(WARNING "NETCDFCXX was NOT found")
	endif()
endif()

# Configure MPI if opted for
if(${WITH_MPI})
	message(STATUS "\nChecking for MPI")
	find_package(MPI QUIET)
	if(MPI_FOUND)
		message(STATUS "MPI was found")
	else()
		message(WARNING "MPI was NOT found")
	endif()
endif()

# Configure GDAL if opted for
if(${WITH_GDAL})
	message(STATUS "\nChecking for GDAL")
	find_package(GDAL REQUIRED QUIET)
	if(GDAL_FOUND)
		message(STATUS "GDAL was found")
	else()
		message(WARNING "GDAL was NOT found")
	endif()
endif()

# Configure CGAL if opted for
if(${WITH_CGAL})
	if(MSVC)
		message(STATUS "\nCGAL is currently disabled on Windows with MSVC as BOOST headers are not compiling")
		#set(CGAL_DISABLE_GMP ON) #GMP is not necessary
		#find_package(CGAL REQUIRED QUIET)
	else()
		message(STATUS "\nChecking for CGAL")
		find_package(CGAL REQUIRED QUIET)
		if(CGAL_FOUND)
			message(STATUS "CGAL was found")
		else()
			message(WARNING "CGAL was NOT found")
		endif()
	endif()
endif()

message(STATUS "== Finished configuring external packages ============================")
message(STATUS "")
