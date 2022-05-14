#!/bin/sh

module load cmake/3.16.2
module load netcdf/4.8.0
module load gdal/3.0.2
module load openmpi/3.0.4

module load intel-compiler
export BUILD_DIR=$PWD/build-intel
export INSTALL_DIR=$PWD/bin/gadi-intel
mkdir $BUILD_DIR
cd $BUILD_DIR
cmake -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DCMAKE_BUILD_TYPE=Release  ..

#module load gcc/11.1.0
#export BUILD_DIR=$PWD/build-gnu
#export INSTALL_DIR=$PWD/bin/gadi-gnu
#mkdir $BUILD_DIR
#cd $BUILD_DIR
#cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_BUILD_TYPE=Release  ..

#module list 

cmake --build . --target all

#cmake --install . --prefix $INSTALL_DIR
#%cmake --install . --prefix /g/data/qi71/apps/utility-programs/bin/gadi/intel
#cmake --install . --prefix /g/data/qi71/apps/utility-programs/bin/gadi/gnu
