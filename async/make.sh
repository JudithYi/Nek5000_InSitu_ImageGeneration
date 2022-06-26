#!/bin/bash
# Add path to Catalyst installation.

#export OSMESA=/path/to/OSMESA

#export LIBDIR=$OSMESA/lib:$LIBDIR
#export LD_LIBRARY_PATH=$OSMESA/lib:$LD_LIBRARY_PATH
#export LD_RUN_PATH=$OSMESA/lib:$LD_RUN_PATH

#export OSMESA_INCLUDE_DIR=$OSMESA/include
#export OSMESA_LIBRARY=$OSMESA/lib

#export PARAVIEW=/path/to/paraview
#export PATH=$PARAVIEW/bin:$PATH
#export LD_LIBRARY_PATH=$PARAVIEW/lib:$LD_LIBRARY_PATH

#export PYTHONPATH=$PARAVIEW/lib/python3.8/site-packages:$PYTHONPATH
#export PYTHONPATH=$PARAVIEW/lib/python3.8/site-packages/paraview/:$PYTHONPATH

export ADIOS2_DIR=/path/to/adios
export PATH=$ADIOS2_DIR/bin:$PATH

export ADIOS2=1
export ADIOS2_INCS=`adios2-config --cxx-flags`
export ADIOS2_LIBS=`adios2-config --cxx-libs`

export CATALYST=0
#export CATALYST_LIBS+=" -L/path/to/paraview -lvtkPVPythonCatalyst-pv5.9 "
#export CATALYST_LIBS+=`paraview-config -l -c Catalyst`
#export CATALYST_INCS=`paraview-config -f -c Catalyst`

export NEK_SOURCE_ROOT="/path/to/Nek5000_Catalyst"

export FC=mpif90
export CC=mpicc
export CXX=mpicxx

${NEK_SOURCE_ROOT}/bin/makenek tgv
./insitu/make.sh
cp insitu/catalystSpace .