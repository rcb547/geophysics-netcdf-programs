#!/bin/sh

module load netcdf/4.7.1
module load gdal/3.0.2
module load openmpi/3.0.4
module load intel-compiler-llvm/2024.0.2
module list 

export BUILD_DIR=$PWD/build-intel
export INSTALL_DIR=$PWD/install-intel
mkdir $BUILD_DIR
cd $BUILD_DIR

# CGAL not currently available on gadi hence -DWITH_CGAL=OFF
cmake -Wno-dev -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx -DWITH_CGAL=OFF ..
cmake --build . --target all
cmake --install . --prefix $INSTALL_DIR
#cmake --install . --prefix /g/data/qi71/apps/utility-programs/bin/gadi/intel

