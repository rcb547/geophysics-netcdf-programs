cmake_minimum_required(VERSION 3.12)

# Set project name
project(geophysics-netcdf-programs)

message(STATUS Build Type = ${CMAKE_BUILD_TYPE})

# Set options
#option(USE_NETCDF "Use NetCDF" ON)
#option(USE_GDAL "Use GDAL" ON)

#message(STATUS USE_NETCDF = ${USE_NETCDF})
#message(STATUS USE_GDAL = ${USE_GDAL})


# Set the C++ standard to C++17
set(CMAKE_CXX_STANDARD 17)

# Configure find package
find_package(PkgConfig REQUIRED)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wno-unknown-pragmas -Wno-unused-variable -Wno-unused-but-set-variable -Wno-sign-compare -Wno-format-security")
SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -O3 -Wno-unknown-pragmas -Wno-unused-variable -Wno-unused-but-set-variable -Wno-sign-compare -Wno-format-security")
if(CMAKE_COMPILER_IS_GNUCC)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-result -Wno-date-time")
	SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -Wno-unused-result -Wno-date-time")
	#On GCC, even with -Wno-date-time, still get warings of the form: warning: macro "__DATE__" might prevent reproducible builds [-Wdate-time]
endif()

# Adding filesystem library
link_libraries(-lstdc++fs)

# Configure GDAL
if(${USE_GDAL})
	pkg_search_module(GDAL gdal IMPORTED_TARGET) 
	if(GDAL_FOUND)
		message(STATUS "GDAL ${GDAL_VERSION} was found")
	endif()
endif()

# Configure NETCDFCXX 
pkg_search_module(NETCDFCXX REQUIRED netcdf-cxx4 IMPORTED_TARGET) 
if(NETCDFCXX_FOUND)
	message(STATUS "NETCDFCXX ${NETCDFCXX_VERSION} was found")
    add_definitions(-DHAVE_NETCDF)
endif()

# Configure MPI
find_package(MPI)
	if(MPI_FOUND)
	message(STATUS "MPI was found")
	include_directories(SYSTEM ${MPI_INCLUDE_PATH})
	add_definitions(-D_MPI_ENABLED -DOMPI_SKIP_MPICXX)
	link_libraries(${MPI_C_LIBRARIES})
endif()


# Add the include directories
include_directories(src/)
include_directories(submodules/cpp-utils/src/)
include_directories(submodules/csv-parser/single_include/)
include_directories(submodules/geophysics-netcdf/src/)
include_directories(submodules/geophysics-netcdf/submodules/marray/include/andres/)

# Add the cpp-utils library submodule
add_subdirectory(submodules/cpp-utils)

# Add the executables
set(mytarget "aseggdf2netcdf")
add_executable(${mytarget}.exe src/${mytarget}.cpp)
target_include_directories(${mytarget}.exe PRIVATE PkgConfig::NETCDFCXX)
target_link_libraries(${mytarget}.exe PRIVATE PkgConfig::NETCDFCXX)
target_link_libraries(${mytarget}.exe PRIVATE cpp-utils-static)
install(TARGETS ${mytarget}.exe DESTINATION .)

set(mytarget "intrepid2netcdf")
add_executable(${mytarget}.exe src/${mytarget}.cpp)
target_include_directories(${mytarget}.exe PRIVATE PkgConfig::NETCDFCXX)
target_link_libraries(${mytarget}.exe PRIVATE PkgConfig::NETCDFCXX)
target_link_libraries(${mytarget}.exe PRIVATE cpp-utils-static)
install(TARGETS ${mytarget}.exe DESTINATION .)

set(mytarget "geophysicsnc2shape")
add_executable(${mytarget}.exe src/${mytarget}.cpp)
target_include_directories(${mytarget}.exe PRIVATE PkgConfig::NETCDFCXX)
target_link_libraries(${mytarget}.exe PRIVATE PkgConfig::NETCDFCXX)
target_link_libraries(${mytarget}.exe PRIVATE cpp-utils-static)
install(TARGETS ${mytarget}.exe DESTINATION .)

set(mytarget "test_geophysics_netcdf")
add_executable(${mytarget}.exe src/${mytarget}.cpp)
target_include_directories(${mytarget}.exe PRIVATE PkgConfig::NETCDFCXX)
target_link_libraries(${mytarget}.exe PRIVATE PkgConfig::NETCDFCXX)
target_link_libraries(${mytarget}.exe PRIVATE cpp-utils-static)
install(TARGETS ${mytarget}.exe DESTINATION .)
