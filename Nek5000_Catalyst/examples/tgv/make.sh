#!/bin/bash
# Add path to Catalyst installation.

export OSMESA=/opt/hpc/build/install

export LIBDIR=$OSMESA/lib:$LIBDIR
export LD_LIBRARY_PATH=$OSMESA/lib:$LD_LIBRARY_PATH
export LD_RUN_PATH=$OSMESA/lib:$LD_RUN_PATH

export OSMESA_INCLUDE_DIR=$OSMESA/include
export OSMESA_LIBRARY=$OSMESA/lib

export PARAVIEW=/opt/hpc/build/install
export PATH=$PARAVIEW/bin:$PATH
export LD_LIBRARY_PATH=$PARAVIEW/lib:$LD_LIBRARY_PATH

export PYTHONPATH=$PARAVIEW/lib/python3.8/site-packages:$PYTHONPATH
export PYTHONPATH=$PARAVIEW/lib/python3.8/site-packages/paraview/:$PYTHONPATH

export ADIOS2_DIR=/opt/hpc/build/install
export PATH=$ADIOS2_DIR/bin:$PATH

export ADIOS2=0
export ADIOS2_INCS=`adios2-config --cxx-flags`
export ADIOS2_LIBS=`adios2-config --cxx-libs`

export MALLEABLE=1
export MALLEABLE_LIBS=$ADIOS2_LIBS

export CATALYST=0
#export CATALYST_LIBS+=" -L/path/to/paraview -lvtkPVPythonCatalyst-pv5.9 "
export CATALYST_LIBS=`paraview-config -l -c PythonCatalyst`
export CATALYST_INCS=`paraview-config -f -c PythonCatalyst`

export MALLEABLE_LIBS+=$CATALYST_LIBS

export LDFLAGS+="-lm -Wl,--copy-dt-needed-entries -lgfortran -L/opt/hpc/install/ompi/lib -Wl,-rpath -Wl,/opt/hpc/install/ompi/lib -Wl,--enable-new-dtags -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi"
export FFLAGS+=-g
export CFLAGS+=-g
export NEK_SOURCE_ROOT=/opt/hpc/build/Nek5000_InSitu_ImageGeneration-main/Nek5000_Catalyst

export FC=mpif90
export CC=mpicc
export CXX=mpicxx

${NEK_SOURCE_ROOT}/bin/makenek tgv

