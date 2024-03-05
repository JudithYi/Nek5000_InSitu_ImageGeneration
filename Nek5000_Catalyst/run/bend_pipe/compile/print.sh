#!/bin/bash
# Add path to Catalyst installation.

export BASEDIR=/u/yju/mpi_session_gcc12_anaconda3_202105
export HWLOCPATH=/u/yju/mpi_session_gcc12_anaconda3_202105/install/hwloc
export LIBEVENTPATH=/u/yju/mpi_session_gcc12_anaconda3_202105/install/libevent
export LIBFABRICPATH=/u/yju/mpi_session_gcc12_anaconda3_202105/install/libfabric

echo "BASE DIR = ${BASEDIR}"
echo "HWLOC PATH = ${HWLOCPATH}"
echo "LIBEVENT PATH = ${LIBEVENTPATH}"
echo "LIBFABRIC_PATH = ${LIBFABRICPATH}"

echo "Exporting the path to the base directory"
export DYNMPI_BASE=$BASEDIR

echo "Exporting the the path to the hwloc and libevent libraries"
export HWLOC_INSTALL_PATH=$HWLOCPATH
export LIBEVENT_INSTALL_PATH=$LIBEVENTPATH
export LIBFABRIC_INSTALL_PATH=$LIBFABRICPATH

echo "Exporting the ompi, open-pmix and prrte install paths"
export PMIX_ROOT=$DYNMPI_BASE/install/pmix
export PRRTE_ROOT=$DYNMPI_BASE/install/prrte
export OMPI_ROOT=$DYNMPI_BASE/install/ompi


export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PMIX_ROOT/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PRRTE_ROOT/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$OMPI_ROOT/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LIBFABRIC_INSTALL_PATH/lib

export C_INCLUDE_PATH=$C_INCLUDE_PATH:$BASEDIR/build/ompi/ompi/include

export PMIX_TOP_SRCDIR=$DYNMPI_BASE/build/openpmix

echo "Updating PATH"
export PATH="$PATH:$OMPI_ROOT/bin"
export PATH="$PATH:$PRRTE_ROOT/bin"
export PATH="$PATH:$LIBFABRIC_INSTALL_PATH/bin"

echo "Environment variables set up successfully"

export OSMESA=/u/yju/mpi_session_gcc12_anaconda3_202105/install/origin/paraview

export LIBDIR=$OSMESA/lib:$LIBDIR
export LD_LIBRARY_PATH=$OSMESA/lib:$LD_LIBRARY_PATH
export LD_RUN_PATH=$OSMESA/lib:$LD_RUN_PATH

export OSMESA_INCLUDE_DIR=$OSMESA/include
export OSMESA_LIBRARY=$OSMESA/lib

export PARAVIEW=/u/yju/mpi_session_gcc12_anaconda3_202105/install/origin/paraview
export PATH=$PARAVIEW/bin:$PATH
export LD_LIBRARY_PATH=$PARAVIEW/lib:$LD_LIBRARY_PATH

export PYTHONPATH=$PARAVIEW/lib/python3.8/site-packages:$PYTHONPATH
export PYTHONPATH=$PARAVIEW/lib/python3.8/site-packages/paraview/:$PYTHONPATH
#export PYTHONPATH=/u/yju/mpi_session_gcc12/python_install/python3.8/site-packages:$PYTHONPATH
export ADIOS2_DIR=/u/yju/mpi_session_gcc12_anaconda3_202105/install/origin/adios
export PATH=$ADIOS2_DIR/bin:$PATH

export ADIOS2=0
export ADIOS2_INCS=`adios2-config --cxx-flags`
export ADIOS2_LIBS=`adios2-config --cxx-libs`

export MALLEABLE=0
export MALLEABLE_LIBS=$ADIOS2_LIBS

export CATALYST=1
#export CATALYST_LIBS+=" -L/path/to/paraview -lvtkPVPythonCatalyst-pv5.9 "
export CATALYST_LIBS=`paraview-config -l -c PythonCatalyst`
export CATALYST_INCS=`paraview-config -f -c PythonCatalyst`

export MALLEABLE_LIBS+=$CATALYST_LIBS

export LDFLAGS+="-lm -Wl,--copy-dt-needed-entries -lgfortran -L/u/yju/mpi_session_gcc12_anaconda3_202105/install/ompi/lib -Wl,-rpath -Wl,/u/yju/mpi_session_gcc12_anaconda3_202105/install/ompi/lib -Wl,--enable-new-dtags -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi"
export FFLAGS+=-g
export CFLAGS+=-g
export NEK_SOURCE_ROOT=/u/yju/mpi_session_gcc12_anaconda3_202105/Nek5000_Catalyst

export FC=/u/yju/mpi_session_gcc12_anaconda3_202105/install/ompi/bin/mpif90
export CC=/u/yju/mpi_session_gcc12_anaconda3_202105/install/ompi/bin/mpicc
export CXX=/u/yju/mpi_session_gcc12_anaconda3_202105/install/ompi/bin/mpicxx
export USR="frame.o mntrlog_block.o mntrlog.o mntrtmr_block.o mntrtmr.o rprm_block.o rprm.o math_tools.o"
USR+=" io_tools_block.o io_tools.o chkpoint.o chkpt_mstp.o gSyEM.o map2D.o stat.o stat_IO.o io_trunc.o"
#${NEK_SOURCE_ROOT}/bin/makenek pipe
echo CATALYST_LIBS
echo $CATALYST_LIBS
echo CATALYST_INCS
echo $CATALYST_INCS

