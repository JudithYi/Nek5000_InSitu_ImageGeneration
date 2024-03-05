### makefile automatically created by makenek 02/12/2024 16:45:43 ###
BINNAME=nek5000
LIBNAME=libnek5000.a
CASENAME:=pipe
CASEDIR:=/u/yju/perf_paper/mix/Nek5000_InSitu_ImageGeneration/Nek5000_Catalyst/run/bend_pipe/compile
S=/u/yju/perf_paper/mix/Nek5000_InSitu_ImageGeneration/Nek5000_Catalyst/core
SCMT=/u/yju/perf_paper/mix/Nek5000_InSitu_ImageGeneration/Nek5000_Catalyst/core/cmt
OPT_INCDIR=./
OBJDIR=obj

FC=mpif90
CC=mpicc
FFLAGS+=-cpp -fdefault-real-8 -fdefault-double-8 -std=legacy -g
CFLAGS+=-g
LDFLAGS+=-lm -Wl,--copy-dt-needed-entries -lgfortran -L/opt/hpc/install/ompi/lib -Wl,-rpath -Wl,/opt/hpc/install/ompi/lib -Wl,--enable-new-dtags -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -Wl,--allow-multiple-definition

PPPO=
PPS= MPI UNDERSCORE GLOBAL_LONG_LONG TIMER MALLEABLE
PPS2= UNDERSCORE GLOBAL_LONG_LONG TIMER

USR=frame.o mntrlog_block.o mntrlog.o mntrtmr_block.o mntrtmr.o rprm_block.o rprm.o math_tools.o io_tools_block.o io_tools.o chkpoint.o chkpt_mstp.o gSyEM.o map2D.o stat.o stat_IO.o io_trunc.o
USR_LFLAGS= -L/u/yju/perf_paper/mix/Nek5000_InSitu_ImageGeneration/Nek5000_Catalyst/3rd_party/blasLapack -lblasLapack -L/u/yju/perf_paper/mix/Nek5000_InSitu_ImageGeneration/Nek5000_Catalyst/3rd_party/gslib/lib -lgs

MPI=1
CMT=0
DPROCMAP=0
VISIT=0

VISIT_IFLAGS=
VISIT_LFLAGS=
HYPRE_IFLAGS=
GSLIB_IFLAGS:= -I/u/yju/perf_paper/mix/Nek5000_InSitu_ImageGeneration/Nek5000_Catalyst/3rd_party/gslib/include
PARRSB_IFLAGS=

ADIOS2:=0
ADIOS2_LFLAGS=
ADIOS2_IFLAGS:=-DADIOS2_USE_MPI -isystem /raven/u/yju/perf_paper/library/install/include 

MALLEABLE:=1
MALLEABLE_LFLAGS:=-Wl,-rpath,/raven/u/yju/perf_paper/library/install/lib64 /raven/u/yju/perf_paper/library/install/lib64/libadios2_cxx11_mpi.so.2.7.1 /raven/u/yju/perf_paper/library/install/lib64/libadios2_cxx11.so.2.7.1 -Wl,-rpath-link,/raven/u/yju/perf_paper/library/install/lib64 -Wl,--no-as-needed -Wl,--unresolved-symbols=ignore-in-shared-libs -Wl,-rpath -Wl,/mpcdf/soft/SLE_15/packages/skylake/openmpi/gcc_10-10.3.0/4.0.7/lib -Wl,--enable-new-dtags -pthread -Wl,-rpath,/u/yju/perf_paper/library/install/lib64:/u/yju/perf_paper/library/install/lib:/raven/u/system/soft/SLE_15/packages/skylake/openmpi/gcc_10-10.3.0/4.0.7/lib:/mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2021.05/lib /u/yju/perf_paper/library/install/lib64/libvtkPVPythonCatalyst-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkWrappingTools-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkRenderingMatplotlib-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkParallelMPI4Py-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkIOMPIImage-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkIOAsynchronous-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkFiltersPython-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkFiltersProgrammable-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkFiltersPoints-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkFiltersParallelStatistics-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkPVCatalyst-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkIOVPIC-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkvpic-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkIOPLY-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkIOParallel-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkIOInfovis-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtklibxml2-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkIOGeometry-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkFiltersParallelVerdict-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkFiltersVerdict-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkverdict-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkFiltersParallelGeometry-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkIOMovie-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkIOParallelXML-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkViewsContext2D-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkRenderingVolumeAMR-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkRenderingVolumeOpenGL2-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkImagingMath-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkRenderingParallel-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkRenderingLabel-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkRenderingContextOpenGL2-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkViewsCore-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkInteractionWidgets-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkImagingHybrid-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkImagingGeneral-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkRenderingAnnotation-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkImagingColor-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkDomainsChemistryOpenGL2-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkDomainsChemistry-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkRenderingOpenGL2-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkglew-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib/libOSMesa.so /u/yju/perf_paper/library/install/lib64/libvtkRenderingUI-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkChartsCore-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkInfovisCore-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkRenderingContext2D-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkInteractionStyle-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkFiltersParallelFlowPaths-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkFiltersFlowPaths-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkRenderingVolume-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkFiltersParallelDIY2-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkFiltersHyperTree-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkFiltersGeneric-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkFiltersParallelMPI-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkFiltersAMR-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkFiltersParallel-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkFiltersTexture-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkFiltersHybrid-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkImagingSources-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkFiltersExtraction-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkParallelDIY-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkjsoncpp-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkPythonInterpreter-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkRenderingFreeType-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkfreetype-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkRenderingCore-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkCommonColor-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkFiltersGeometry-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkParallelMPI-pv5.9.so.5.9 /raven/u/system/soft/SLE_15/packages/skylake/openmpi/gcc_10-10.3.0/4.0.7/lib/libmpi.so /u/yju/perf_paper/library/install/lib64/libvtkIOXML-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkIOXMLParser-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkexpat-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkIOImage-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkDICOMParser-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkmetaio-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkpng-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtktiff-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkjpeg-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkWrappingPythonCore3.8-pv5.9.so.5.9 /mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2021.05/lib/libpython3.8.so /u/yju/perf_paper/library/install/lib64/libvtkFiltersModeling-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkFiltersSources-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkFiltersGeneral-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkCommonComputationalGeometry-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkFiltersCore-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkParallelCore-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkIOLegacy-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkIOCore-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtklzma-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtklz4-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkdoubleconversion-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkzlib-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkFiltersStatistics-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkImagingFourier-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkImagingCore-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkCommonExecutionModel-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkCommonDataModel-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkCommonTransforms-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkpugixml-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkCommonMisc-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkCommonMath-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkCommonSystem-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkCommonCore-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtkloguru-pv5.9.so.5.9 /u/yju/perf_paper/library/install/lib64/libvtksys-pv5.9.so.5.9 -ldl -lpthread -Wl,-rpath-link,/u/yju/perf_paper/library/install/lib64 -lstdc++

CATALYST:=0
CATALYST_LFLAGS=
CATALYST_IFLAGS:=-DDIY_NO_THREADS -DMPICH_SKIP_MPICXX -DOMPI_SKIP_MPICXX -D_MPICC_H -isystem /u/yju/perf_paper/library/install/include/paraview-5.9 -isystem /mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2021.05/include/python3.8 -isystem /u/yju/perf_paper/library/install/include -isystem /u/yju/perf_paper/library/install/include/paraview-5.9/vtkmpi4py -isystem /u/yju/perf_paper/library/install/include/paraview-5.9/vtklibxml2/include -isystem /raven/u/system/soft/SLE_15/packages/skylake/openmpi/gcc_10-10.3.0/4.0.7/include -isystem /raven/u/system/soft/SLE_15/packages/skylake/openmpi/gcc_10-10.3.0/4.0.7/include/openmpi -isystem /raven/u/system/soft/SLE_15/packages/skylake/openmpi/gcc_10-10.3.0/4.0.7/include/openmpi/opal/mca/hwloc/hwloc201/hwloc/include -isystem /raven/u/system/soft/SLE_15/packages/skylake/openmpi/gcc_10-10.3.0/4.0.7/include/openmpi/opal/mca/event/libevent2022/libevent -isystem /raven/u/system/soft/SLE_15/packages/skylake/openmpi/gcc_10-10.3.0/4.0.7/include/openmpi/opal/mca/event/libevent2022/libevent/include -pthread


################################################################################

PPS_F = $(patsubst %,$(PPPO)-D%,$(PPS))
PPS_C = $(patsubst %,-D%,$(PPS))
PPS_C2 = $(patsubst %,-D%,$(PPS2))

# CORE##########################################################################
CORE = drive1_insitu.o drive2.o \
plan5.o \
plan4.o bdry.o coef.o conduct.o connect1.o connect2.o \
dssum.o edgec.o eigsolv.o gauss.o genxyz.o navier1.o makeq.o \
navier0.o navier2.o navier3.o navier4.o prepost.o speclib.o \
map2.o mvmesh.o ic.o gfldr.o ssolv.o planx.o math.o mxm_wrapper.o \
hmholtz.o subs1.o subs2.o \
gmres.o hsmg.o convect.o convect2.o induct.o perturb.o \
navier5.o navier6.o navier7.o navier8.o fast3d.o fasts.o calcz.o \
byte.o chelpers.o byte_mpi.o postpro.o interp.o \
cvode_driver.o multimesh.o \
vprops.o makeq_aux.o \
papi.o nek_in_situ.o \
reader_rea.o reader_par.o reader_re2.o \
finiparser.o iniparser.o dictionary.o \
hpf.o \
fcrs.o crs_xxt.o crs_amg.o \
fem_amg_preco.o crs_hypre.o \
partitioner.o
################################################################################
# MXM 
MXM=mxm_std.o
################################################################################

COMM_MPI := comm_mpi.o
ifeq ($(MPI),0)
  COMM_MPI := ${COMM_MPI} mpi_dummy.o
endif

# CMT append *.o files to CORE
ifeq ($(CMT),0)
  CMTCORE := 
else
  CMTCORE := drive1_cmt.o drive2_cmt.o driver3_cmt.o \
   face.o ausm.o  state1.o intpdiff.o diffusive_cmt.o eqnsolver_cmt.o surface_fluxes.o \
   bc.o inflow_bc.o outflow_bc.o wall_bc.o step.o \
   filters_cmt.o diagnostics.o artvisc.o
endif

DUMMY:= $(shell cp $S/PARALLEL.default $S/PARALLEL)
ifeq ($(DPROCMAP),1)
	CORE := ${CORE} dprocmap.o
 	DUMMY:= $(shell cp $S/PARALLEL.dprocmap $S/PARALLEL)
endif

ifneq ($(VISIT),0)
  VISITNEK_INCLUDES:=-I$S/3rd_party
  VISITO := visit.o 
  USR_LFLAGS := ${USR_LFLAGS} $(VISIT_LFLAGS)
endif

# ADIOS2 #####################################################################
ifneq ($(ADIOS2),0)
  ADIOS2O := adios2.o nek_adios2.o 
  USR_LFLAGS := ${USR_LFLAGS} $(ADIOS2_LFLAGS)
endif

# CATALYST #####################################################################
ifneq ($(CATALYST),0)
  CATALYSTO := catalyst.o nek_catalyst.o 
  USR_LFLAGS := ${USR_LFLAGS} $(CATALYST_LFLAGS)
endif

# MALLEABLE #####################################################################
ifneq ($(MALLEABLE),0)
  MALLEABLEO := adios2.o nek_adios2.o adios_catalyst_malleable_insitu.o myCPPythonAdaptorAPI.o nek_catalyst_async.o
  USR_LFLAGS := ${USR_LFLAGS} $(MALLEABLE_LFLAGS)
endif

ifeq ($(MPI),0)
	DUMMY:= $(shell cp $S/mpi_dummy.h $S/mpif.h) 
else
	DUMMY:= $(shell rm -rf $S/mpif.h) 
endif

TMP1 = $(CORE) $(MXM) $(USR) $(COMM_MPI) $(VISITO) $(ADIOS2O) $(CATALYSTO) $(MALLEABLEO)

ifeq ($(CMT),1)
TMP1 += $(CMTCORE)
OPT_INCDIR += -I $(SCMT)
endif

OPT_INCDIR += -I $S/experimental

NOBJS = $(patsubst %,$(OBJDIR)/%,$(TMP1))

L0 = $(G) -O0
L2 = $(G) -O2
L3 = $(G) -O3

FL0   = $(L0) $(FFLAGS) $(PPS_F) -I$(CASEDIR) -I$S -I$(OPT_INCDIR) 
FL2   = $(L2) $(FFLAGS) $(PPS_F) -I$(CASEDIR) -I$S -I$(OPT_INCDIR)
FL3   = $(L3) $(FFLAGS) $(PPS_F) -I$(CASEDIR) -I$S -I$(OPT_INCDIR)

cFL0   = $(L0) $(CFLAGS) $(PPS_C) -I$S -I$(OPT_INCDIR)
cFL2   = $(L2) $(CFLAGS) $(PPS_C) -I$S -I$(OPT_INCDIR)
cFL22   = $(L2) $(CFLAGS) $(PPS_C2) -I$S -I$(OPT_INCDIR)
cFL3   = $(L3) $(CFLAGS) $(PPS_C) -I$S -I$(OPT_INCDIR)

################################################################################
.PHONY: all clean

all : nek5000

objdir: 
	@mkdir $(OBJDIR) 2>/dev/null; cat /dev/null 

usrfile:
	@env CASENAME=$(CASENAME) PPS="$(PPS)" $S/mkuserfile
	$(FC) $(FL2) $(VISITNEK_INCLUDES) -c $(CASEDIR)/${CASENAME}.f

nek5000: lib usrfile
	$(CXX) $(cFL2) $(ADIOS2_IFLAGS) $(CATALYST_IFLAGS) -o ${BINNAME} $S/drive.cpp $(CASEDIR)/${CASENAME}.o ${LIBNAME} $(USR_LFLAGS) $(LDFLAGS)
#	$(FC) $(FL2) -o ${BINNAME} $S/drive.f $(CASEDIR)/${CASENAME}.o ${LIBNAME} $(USR_LFLAGS) $(LDFLAGS)
	@if test -f ${BINNAME}; then \
        echo ""; \
        echo "$(NEK_WARN)"; \
	echo "#############################################################"; \
	echo "#                  Compilation successful!                  #"; \
	echo "#############################################################"; \
	size ${BINNAME}; \
	echo ""; \
	else \
	echo -e "ERROR: Cannot find ${BINNAME}!"; \
	echo -e ""; \
	fi
	@rm -f ${CASENAME}.o ${LIBNAME} drive.o *.i

lib:	objdir $(NOBJS)
	@printf "building ${LIBNAME} ... "
	@$(AR) cru ./${LIBNAME} $(NOBJS)
	@ranlib ${LIBNAME}
	@if test -f ${LIBNAME}; then \
	printf "done\n"; \
	else \
	printf "failed!\n"; \
	echo ""; \
	fi
	@rm -f *.i

clean::
	@rm -rf ${LIBNAME}  ${CASENAME}.o ${CASENAME}.f obj ${BINNAME}

$(NOBJS) : ./SIZE

# NEK CORE     ##################################################################
$(OBJDIR)/drive1.o	:$S/drive1.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/drive1_insitu.o	:$S/drive1_insitu.f;	$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/drive2.o	:$S/drive2.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/prepost.o	:$S/prepost.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/postpro.o	:$S/postpro.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/gfldr.o	:$S/gfldr.f $S/GFLDR;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/connect1.o	:$S/connect1.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/connect2.o	:$S/connect2.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/edgec.o	:$S/edgec.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/genxyz.o	:$S/genxyz.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/subs1.o	:$S/subs1.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/subs2.o	:$S/subs2.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/map2.o	:$S/map2.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/dprocmap.o    :$S/dprocmap.f $S/DPROCMAP;     $(FC) -c $(FL2) $< -o $@
$(OBJDIR)/interp.o      :$S/interp.f;                   $(FC) -c $(FL2) $< -o $@
$(OBJDIR)/mvmesh.o	:$S/mvmesh.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/induct.o	:$S/induct.f;	      		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/convect.o	:$S/convect.f;	      		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/convect2.o	:$S/convect2.f;	      		$(FC) -c $(FL0) $< -o $@
$(OBJDIR)/perturb.o	:$S/perturb.f;	      		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/hsmg.o	:$S/hsmg.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/gmres.o	:$S/gmres.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/papi.o	:$S/papi.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/bdry.o	:$S/bdry.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/ic.o		:$S/ic.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/conduct.o	:$S/conduct.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/navier0.o	:$S/navier0.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/navier2.o	:$S/navier2.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/navier3.o	:$S/navier3.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/navier4.o	:$S/navier4.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/navier5.o	:$S/navier5.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/navier6.o	:$S/navier6.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/navier7.o	:$S/navier7.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/navier8.o	:$S/navier8.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/speclib.o	:$S/speclib.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/calcz.o	:$S/calcz.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/fast3d.o	:$S/fast3d.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/dssum.o	:$S/dssum.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/eigsolv.o	:$S/eigsolv.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/gauss.o	:$S/gauss.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/planx.o	:$S/planx.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/ssolv.o	:$S/ssolv.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/coef.o	:$S/coef.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/plan4.o	:$S/plan4.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/plan5.o	:$S/plan5.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/vprops.o	:$S/vprops.f;	   		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/cvode_driver.o :$S/cvode_driver.f $S/CVODE;	$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/makeq.o      	:$S/makeq.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/makeq_aux.o  	:$S/makeq_aux.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/hmholtz.o	:$S/hmholtz.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/navier1.o	:$S/navier1.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/fasts.o	:$S/fasts.f;		        $(FC) -c $(FL2) $< -o $@
$(OBJDIR)/comm_mpi.o	:$S/comm_mpi.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/mpi_dummy.o	:$S/mpi_dummy.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/byte_mpi.o	:$S/byte_mpi.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/math.o	:$S/math.f;			$(FC) -c $(FL3) $< -o $@
$(OBJDIR)/multimesh.o	:$S/multimesh.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/lb_setqvol.o	:$S/lb_setqvol.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/reader_rea.o	:$S/reader_rea.f;	 	$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/reader_par.o	:$S/reader_par.f $S/PARDICT;	$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/reader_re2.o	:$S/reader_re2.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/hpf.o	       	:$S/hpf.f;         		$(FC) -c $(FL2) $< -o $@

# CMT CORE     ##################################################################
$(OBJDIR)/step.o	:$(SCMT)/step.f;        	$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/drive1_cmt.o	:$(SCMT)/drive1_cmt.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/drive2_cmt.o	:$(SCMT)/drive2_cmt.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/driver3_cmt.o	:$(SCMT)/driver3_cmt.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/gradient_solver.o	:$(SCMT)/gradient_solver.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/eqnsolver_cmt.o	:$(SCMT)/eqnsolver_cmt.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/ausm.o	:$(SCMT)/ausm.f;   		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/artvisc.o	:$(SCMT)/artvisc.f;   		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/face.o	:$(SCMT)/face.f;  		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/bc.o		:$(SCMT)/bc.f;			$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/surface_fluxes.o	:$(SCMT)/surface_fluxes.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/diffusive_cmt.o	:$(SCMT)/diffusive_cmt.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/state1.o	:$(SCMT)/state1.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/intpdiff.o	:$(SCMT)/intpdiff.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/inflow_bc.o	:$(SCMT)/inflow_bc.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/outflow_bc.o	:$(SCMT)/outflow_bc.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/wall_bc.o	:$(SCMT)/wall_bc.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/filters_cmt.o	:$(SCMT)/filters_cmt.f;		$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/diagnostics.o	:$(SCMT)/diagnostics.f;		$(FC) -c $(FL2) $< -o $@

# MXM       ############################################################################
$(OBJDIR)/mxm_wrapper.o	  :$S/mxm_wrapper.f;		$(FC) -c $(FL2) $< -o $@ 
$(OBJDIR)/mxm_std.o	  :$S/mxm_std.f;		$(FC) -c $(FL3) $< -o $@
$(OBJDIR)/mxm_bgq.o	  :$S/mxm_bgq.f;		$(FC) -c $(FL3) $< -o $@

# C Files ##################################################################################
$(OBJDIR)/byte.o             :$S/byte.c;              $(CC) -c $(cFL2) $< -o $@
$(OBJDIR)/chelpers.o         :$S/chelpers.c;          $(CC) -c $(cFL2) $< -o $@
$(OBJDIR)/fcrs.o             :$S/fcrs.c;              $(CC) -c $(cFL2) $(GSLIB_IFLAGS) $(HYPRE_IFLAGS) $< -o $@
$(OBJDIR)/crs_xxt.o          :$S/crs_xxt.c;           $(CC) -c $(cFL2) $(GSLIB_IFLAGS) $< -o $@
$(OBJDIR)/crs_amg.o          :$S/crs_amg.c;           $(CC) -c $(cFL2) $(GSLIB_IFLAGS) $< -o $@
$(OBJDIR)/fem_amg_preco.o    :$S/experimental/fem_amg_preco.c; $(CC) -c $(cFL2) $(GSLIB_IFLAGS) $(HYPRE_IFLAGS) $< -o $@
$(OBJDIR)/crs_hypre.o        :$S/experimental/crs_hypre.c; $(CC) -c $(cFL2) $(GSLIB_IFLAGS) $(HYPRE_IFLAGS) $< -o $@
$(OBJDIR)/partitioner.o      :$S/partitioner.c;     $(CC) -c $(cFL2) $(GSLIB_IFLAGS) $(PARRSB_IFLAGS) $< -o $@

# 3rd party #######################################################################################
$(OBJDIR)/visit.o        :$S/3rd_party/visit.f;		$(FC) -c $(VISIT_IFLAGS) $(FL2) $< -o $@
$(OBJDIR)/adios2.o     :$S/3rd_party/adios2.f;	$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/nek_adios2.o :$S/3rd_party/nek_adios2.cpp;$(CXX)  -c $(cFL22) $(ADIOS2_IFLAGS) $< -o $@
$(OBJDIR)/catalyst.o     :$S/3rd_party/catalyst.f;	$(FC) -c $(FL2) $< -o $@
$(OBJDIR)/myCPPythonAdaptorAPI.o :$S/3rd_party/myCPPythonAdaptorAPI.cpp;$(CXX)  -c $(cFL22) $(CATALYST_IFLAGS) $< -o $@
$(OBJDIR)/nek_catalyst_async.o :$S/3rd_party/nek_catalyst_async.cpp;$(CXX)  -c $(cFL22) $(ADIOS2_IFLAGS) $(CATALYST_IFLAGS) $< -o $@
$(OBJDIR)/adios_catalyst_malleable_insitu.o :$S/3rd_party/adios_catalyst_malleable_insitu.cpp;$(CXX)  -c $(cFL22) $(ADIOS2_IFLAGS) $(CATALYST_IFLAGS) $< -o $@
$(OBJDIR)/nek_catalyst.o :$S/3rd_party/nek_catalyst.cxx;$(CXX)  -c $(cFL2) $(CATALYST_IFLAGS) $< -o $@
$(OBJDIR)/nek_in_situ.o  :$S/3rd_party/nek_in_situ.f;	$(FC) -c $(VISIT_IFLAGS) $(FL2) $< -o $@
$(OBJDIR)/finiparser.o   :$S/3rd_party/finiparser.c;    $(CC)  -c $(cFL2) $(JL) $< -o $@
$(OBJDIR)/iniparser.o    :$S/3rd_party/iniparser.c;     $(CC)  -c $(cFL2) $(JL) $< -o $@
$(OBJDIR)/dictionary.o   :$S/3rd_party/dictionary.c;    $(CC)  -c $(cFL2) $(JL) $< -o $@
include makefile_usr.inc
