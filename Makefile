# Makefile for various platforms
# Execute using Build csh-script only!
# Used together with Perl scripts in SRC/SCRIPT 
# (C) 2005 Marat Khairoutdinov
#------------------------------------------------------------------
# uncomment to disable timers:
#
#NOTIMERS=-DDISABLE_TIMERS
#-----------------------------------------------------------------

SAM = SAM_$(RAD_DIR)_$(MICRO_DIR)_w_aero

# Determine platform 
PLATFORM := $(shell uname -s)

#----------------------------------
#HOSTNAME = $(shell hostname)
#ifeq ($(HOSTNAME),rotor)
  TUV_SRC_DIR = /nas/data_cellar/rotor/data/p2033/SAM-ASP_tuv/SRC_TUV/TUV

  ASP_SRC_DIR = /nas/data_cellar/rotor/data/p2033/SAM-ASP_tuv/SRC_TUV/ASP

  ASP_LIBRARY_DIR = $(ASP_SRC_DIR)

  ##rotor##LIB_MPI = /rotor/data/p1808/CMAQv5.0.1/lib/x86_64/ifort/mpich/lib -lmpi -lrdmacm -libumad -lrt -lpthread -libverbs -ldl -lxmpi
  LIB_MPI = /usr/lib64/mpich-3.2/lib -lmpi -lrdmacm -libumad -lrt -lpthread -libverbs -ldl -lmpifort
  ##rotor##INC_MPI =  /rotor/data/p1808/CMAQv5.0.1/lib/x86_64/ifort/mpich/include
  INC_MPI = /usr/include/mpich-3.2-x86_64
  ##rotor##INC_NETCDF = /rotor/data/p1808/CMAQv5.0.1/lib/x86_64/ifort/netcdf/include
  INC_NETCDF = /usr/include

  ##rotor##LIB_NETCDF = /rotor/data/p1808/CMAQv5.0.1/lib/x86_64/ifort/netcdf/lib -lnetcdff -lnetcdf -lnetcdf
  LIB_NETCDF = /usr/lib64 -lnetcdff -lnetcdf -lnetcdf

  ##rotor##HDF5_HL =  /usr/local/hdf5/1.8.6/intel/lib -lhdf5_hl -lhdf5
  HDF5_HL = /usr/lib64 -lhdf5_hl -lhdf5

  ##rotor##Z =  /usr/local/zlib/1.2.5/intel/lib -lz -lm -lhdf5_hl -lhdf5 -lz
  Z =  /usr/lib64 -lz -lm -lhdf5_hl -lhdf5 -lz

  FF77 = ifort -c -fixed -extend_source
  FF90 = ifort -c -free
  FFLAGS = -I$(INC_NETCDF) -I$(INC_MPI) -g -O2 -pad -I$(ASP_LIBRARY_DIR)/ModuleFiles 
#  FFLAGS = -I$(INC_NETCDF) -I$(INC_MPI) -g -O2 -pad -fpe0 -g -traceback -debug extended -I$(ASP_LIBRARY_DIR)/ModuleFiles #AD DEBUG

  LD = ifort
  LDFLAGS = -L$(LIB_MPI) -L$(LIB_NETCDF) -L$(HDF5_HL) -L$(Z)
  LDFLAGS += -L${ASP_SRC_DIR} -lAspLibrary -L${TUV_SRC_DIR} -lTuvLibrary
  CC = icc -c -DLINUX
#endif 


#compute the search path
dirs := . $(shell cat Filepath)
VPATH    := $(foreach dir,$(dirs),$(wildcard $(dir))) 

.SUFFIXES:
.SUFFIXES: .f .f90 .c .o



all: $(SAM_DIR)/$(SAM)
#all:
#	make Srcfiles
#	make Depends
#	make $(SAM_DIR)/$(SAM)

SOURCES   := $(shell cat Srcfiles)

Depends: Srcfiles Filepath
	perl $(SAM_SRC)/SCRIPT/mkDepends Filepath Srcfiles > $@

Srcfiles: Filepath
	perl $(SAM_SRC)/SCRIPT/mkSrcfiles > $@

OBJS      := $(addsuffix .o, $(basename $(SOURCES))) 

$(SAM_DIR)/$(SAM): $(OBJS)
	$(LD) -o $@ $(OBJS) $(LDFLAGS)


.f90.o:
	${FF90}  ${FFLAGS} $<
.f.o:
	${FF77}  ${FFLAGS} $<
.c.o:
	${CC}  ${CFLAGS} -I$(SAM_SRC)/TIMING $(NOTIMERS) $<

include Depends

clean: 
	rm $(SAM_DIR)/OBJ/*


