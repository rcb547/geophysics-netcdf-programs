#!/bin/sh

module load netcdf/4.7.1
module load gdal/3.0.2
module load openmpi/3.0.4
module load gcc/13.2.0
module list 

export BUILD_DIR=$PWD/build-gnu
export INSTALL_DIR=$PWD/install-gnu
mkdir $BUILD_DIR
cd $BUILD_DIR

# CGAL not currently available on gadi hence -DWITH_CGAL=OFF
cmake -Wno-dev -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DWITH_CGAL=OFF ..
cmake --build . --target all
cmake --install . --prefix $INSTALL_DIR
#cmake --install . --prefix /g/data/qi71/apps/utility-programs/bin/gadi/gnu
