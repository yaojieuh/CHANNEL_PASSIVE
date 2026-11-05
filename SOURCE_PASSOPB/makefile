
MACHINE = ZHUHAI
#MACHINE = TTU_impi
LDFLAGS =
SOURCES = module_timer.f90					\
	module_global_variable.f90				\
	module_labs.f90						\
	module_runge_kutta.f90					\
	module_ppdfftw.f90					\
	print_program_info.f90					\
	check_remaining_time.f90				\
	labs_fill_full.cpp					\
	read_input_file.f90					\
	wavespace.f90						\
	dft_x_realspace.f90					\
	x_to_y.f90						\
	y_to_x.f90						\
	wave00.f90						\
	initial_step.f90					\
	get_statistics.f90					\
	get_pressure.f90					\
	write_energy_spectra.f90				\
	data_store.f90						\
	rescale_data.f90					\
	solve_v_f.f90						\
	solve_g.f90						\
	solve_t.f90						\
	solve_u_w.f90						\
	hvhg.f90						\
	time_march.f90						\
	main.f90

EXEC    = POONGBACK_Channel_$(MACHINE)

##############################################################
## building object list ##

OTHERLIBS=
OTHERINC=

OBJECTS = $(addsuffix .o,$(basename $(SOURCES)))
########################################################################
ifeq ($(MACHINE),ZHUHAI)
COMPILE = mpif90 -O3 -fopenmp
CC= mpicc -O3 -fopenmp -fcommon
#FC= mpif90 -O3 -fopenmp -fallow-argument-mismatch -fcommon
FC= mpif90 -O3 -fopenmp  -fcommon
GSL_LIB=/home/yaoj/SOURCE_FILE/ESIO/gsl-1.16-build/lib
GSL_DIR=/home/yaoj/SOURCE_FILE/ESIO/gsl-1.16-build
#GSL_LIB=/opt/apps/gcc13/gsl/2.8/lib
#GSL_DIR=/opt/apps/gcc13/gsl/2.8/
FFTW_LIB=/home/yaoj/SOURCE_FILE/ESIO/fftw-3.3.7-build/lib
FFTW_INC=/home/yaoj/SOURCE_FILE/ESIO/fftw-3.3.7-build/include
HDF_DIR=/home/yaoj/SOURCE_FILE//ESIO/hdf5-1.8.18-build
ESIO_LIB=/home/yaoj/SOURCE_FILE//ESIO/ESIO-release-0.1.9-build/lib
ESIO_INC=/home/yaoj/SOURCE_FILE/ESIO/ESIO-release-0.1.9-build/include

GSL_LIBRARY=-L$(GSL_LIB) -lgsl -lgslcblas -lstdc++ -lm
GSL_INCLUDE = -I$(GSL_DIR)/include
FFTW3_LIB=-L$(FFTW_LIB) -lfftw3_mpi -lfftw3 -lm
FFTW3_INC=-I$(FFTW_INC)
ESIO_LIBRARY = -L$(ESIO_LIB) -lesiof -lm
ESIO_INCLUDE = -I$(ESIO_LIB)
HDF_LIB=-L$(HDF_DIR)/lib -lhdf5_hl -lhdf5 -lm -lz
HDF_INC=-I$(HDF_DIR)/include
endif
########################################################################
########################################################################
ifeq ($(MACHINE),TACC)
COMPILE = mpif90 -O3 -fopenmp
CC= mpicc -O3 -fopenmp -fcommon
FC= mpif90 -O3 -fopenmp -fallow-argument-mismatch -fcommon
GSL_LIB=/home1/05779/fazlehus/ESIO/gsl-2.7.1-build/lib
GSL_DIR=/home1/05779/fazlehus/ESIO/gsl-2.7.1-build
FFTW_LIB=/opt/apps/gcc13/impi21/fftw3/3.3.10/lib
FFTW_INC=/opt/apps/gcc13/impi21/fftw3/3.3.10//include
HDF_DIR=/home1/05779/fazlehus/ESIO/hdf5-1.8.18-build
ESIO_LIB=/home1/05779/fazlehus/ESIO/ESIO-release-0.1.9-build/lib
ESIO_INC=/home1/05779/fazlehus/ESIO/ESIO-release-0.1.9-build/include

GSL_LIBRARY=-L$(GSL_LIB) -lgsl -lgslcblas -lstdc++ -lm
GSL_INCLUDE = -I$(GSL_DIR)/include
FFTW3_LIB=-L$(FFTW_LIB) -lfftw3_mpi -lfftw3 -lm
FFTW3_INC=-I$(FFTW_INC)
ESIO_LIBRARY = -L$(ESIO_LIB) -lesiof -lm
ESIO_INCLUDE = -I$(ESIO_LIB)
HDF_LIB=-L$(HDF_DIR)/lib -lhdf5_hl -lhdf5 -lm -lz
HDF_INC=-I$(HDF_DIR)/include
endif
########################################################################
ifeq ($(MACHINE),TTU_impi)
COMPILE = mpif90 -O3 -fopenmp
CC= mpicc -O3 -fopenmp
FC= mpif90 -O3 -fopenmp
GSL_LIB=/lustre/work/jieyao/ESIO-ii/gsl-1.16-build/lib
GSL_DIR=/lustre/work/jieyao/ESIO-ii/gsl-1.16-build
FFTW_LIB= /lustre/work/jieyao/ESIO-ii/fftw-3.3.7-build/lib
FFTW_INC=/lustre/work/jieyao/ESIO-ii/fftw-3.3.7-build/include
#HDF_DIR=/lustre/work/jieyao/esio-019/hdf5-1.10.0-build
#HDF_DIR=/lustre/work/apps/hdf5-1.8.18-parallel
HDF_DIR=/lustre/work/jieyao/ESIO-ii/hdf5-1.8.18-build/
ESIO_LIB=/lustre/work/jieyao/ESIO-ii/esio-0.1.9-build/lib
ESIO_INC=/lustre/work/jieyao/ESIO-ii/esio-0.1.9-build/include
GSL_LIBRARY=-L$(GSL_LIB) -lgsl -lgslcblas -lstdc++ -lm
GSL_INCLUDE = -I$(GSL_DIR)/include
FFTW3_LIB=-L$(FFTW_LIB) -lfftw3_mpi -lfftw3 -lm
FFTW3_INC=-I$(FFTW_INC)
ESIO_LIBRARY = -L$(ESIO_LIB) -lesiof -lm
ESIO_INCLUDE = -I$(ESIO_LIB)
HDF_LIB=-L$(HDF_DIR)/lib -lhdf5_hl -lhdf5 -lm -lz
HDF_INC=-I$(HDF_DIR)/include
endif
########################################################################



########################################################################
ifeq ($(MACHINE),opuntia)
COMPILE = mpif90 -O3 -fopenmp
CC= mpicc -O3 -fopenmp
FC= mpif90 -O3 -fopenmp
GSL_LIB= /project/kouri/ESIO/gsl-1.16-build/lib
GSL_DIR=/project/kouri/ESIO/gsl-1.16-build
#FFTW_LIB= /project/cacds/apps/fftw/3.3.4/lib/
#FFTW_INC=/project/cacds/apps/fftw/3.3.4/include
FFTW_LIB=/project/kouri/ESIO/fftw-3.3.7-build/lib
FFTW_INC=/project/kouri/ESIO/fftw-3.3.7-build/include

HDF_DIR= /project/kouri/ESIO/hdf5-1.8.18-build
ESIO_LIB=/project/kouri/ESIO/ESIO-release-0.1.9-build/lib
ESIO_INC=/project/kouri/ESIO/ESIO-release-0.1.9-build/lib
GSL_LIBRARY=-L$(GSL_LIB) -lgsl -lgslcblas -lstdc++ -lm
GSL_INCLUDE = -I$(GSL_DIR)/include
FFTW3_LIB=-L$(FFTW_LIB) -lfftw3_mpi -lfftw3 -lm
FFTW3_INC=-I$(FFTW_INC)
ESIO_LIBRARY = -L$(ESIO_LIB) -lesiof -lm
ESIO_INCLUDE = -I$(ESIO_INC)	
HDF_LIB=-L$(HDF_DIR)/lib -lhdf5_hl -lhdf5 -lm -lz
HDF_INC=-I$(HDF_DIR)/include
endif
########################################################################

########################################################################
ifeq ($(MACHINE),nocona)
COMPILE = mpifort -O3 -fopenmp
CC= mpicc -O3 -fopenmp -fcommon
FC= mpif90 -O3 -fopenmp -fallow-argument-mismatch -fcommon
GSL_LIB=/home/jieyao/CODE_TEST/CHANNEL/gsl-1.16-build/lib
GSL_DIR=/home/jieyao/CODE_TEST/CHANNEL/gsl-1.16-build
#FFTW_LIB=/opt/apps/nfs/spack/opt/spack/linux-centos8-zen2/gcc-10.1.0/fftw-3.3.8-ixu7d5257fwbh2fptwqgne4qso5rva5v/lib
#FFTW_INC=/opt/apps/nfs/spack/opt/spack/linux-centos8-zen2/gcc-10.1.0/fftw-3.3.8-ixu7d5257fwbh2fptwqgne4qso5rva5v/include
FFTW_LIB=/home/jieyao/CODE_TEST/CHANNEL/fftw-3.3.7-build/lib
FFTW_INC=/home/jieyao/CODE_TEST/CHANNEL/fftw-3.3.7-build//include

#HDF_DIR=/opt/apps/nfs/spack/opt/spack/linux-centos8-zen2/gcc-10.1.0/hdf5-1.10.6-bfqfwx3iaupxuyeqzx7lshc424gfeang
HDF_DIR=/home/jieyao/CODE_TEST/CHANNEL/hdf5-1.8.18-build
ESIO_LIB=/home/jieyao/CODE_TEST/CHANNEL/esio-0.1.9-build/lib
ESIO_INC=/home/jieyao/CODE_TEST/CHANNEL/esio-0.1.9-build/include
GSL_LIBRARY=-L$(GSL_LIB) -lgsl -lgslcblas -lstdc++ -lm
GSL_INCLUDE = -I$(GSL_DIR)/include
FFTW3_LIB=-L$(FFTW_LIB) -lfftw3_mpi -lfftw3 -lm
FFTW3_INC=-I$(FFTW_INC)
ESIO_LIBRARY = -L$(ESIO_LIB) -lesiof -lm
ESIO_INCLUDE = -I$(ESIO_LIB)
HDF_LIB=-L$(HDF_DIR)/lib -lhdf5_hl -lhdf5 -lm -lz
HDF_INC=-I$(HDF_DIR)/include
endif
########################################################################

########################################################################
ifeq ($(MACHINE),nocona2)
COMPILE = mpifort -O3 -fopenmp
CC= mpicc -O3 -fopenmp -fcommon
FC= mpif90 -O3 -fopenmp -fallow-argument-mismatch -fcommon
GSL_LIB=/home/jieyao/CODE_TEST/CHANNEL/gsl-1.16-build/lib
GSL_DIR=/home/jieyao/CODE_TEST/CHANNEL/gsl-1.16-build
#FFTW_LIB=/opt/apps/nfs/spack/opt/spack/linux-centos8-zen2/gcc-10.1.0/fftw-3.3.8-dja5okiklo4hh2pifkiaonqkv2lmcjda/lib
#FFTW_INC=/opt/apps/nfs/spack/opt/spack/linux-centos8-zen2/gcc-10.1.0/fftw-3.3.8-dja5okiklo4hh2pifkiaonqkv2lmcjda/include
FFTW_LIB=/lustre/work/thomasbr/spack-yao/opt/spack/linux-centos8-zen2/gcc-10.1.0/fftw-3.3.8-kxrcyt6kfhcnepiotsg6pxraeoecxays/lib
FFTW_INC=/lustre/work/thomasbr/spack-yao/opt/spack/linux-centos8-zen2/gcc-10.1.0/fftw-3.3.8-kxrcyt6kfhcnepiotsg6pxraeoecxays/include

#HDF_DIR=/opt/apps/nfs/spack/opt/spack/linux-centos8-zen2/gcc-10.1.0/hdf5-1.10.6-bfqfwx3iaupxuyeqzx7lshc424gfeang
HDF_DIR=/home/jieyao/CODE_TEST/CHANNEL/hdf5-1.8.18-build
ESIO_LIB=/home/jieyao/CODE_TEST/CHANNEL/esio-0.1.9-build/lib
ESIO_INC=/home/jieyao/CODE_TEST/CHANNEL/esio-0.1.9-build/include
GSL_LIBRARY=-L$(GSL_LIB) -lgsl -lgslcblas -lstdc++ -lm
GSL_INCLUDE = -I$(GSL_DIR)/include
FFTW3_LIB=-L$(FFTW_LIB) -lfftw3_mpi -lfftw3 -lm
FFTW3_INC=-I$(FFTW_INC)
ESIO_LIBRARY = -L$(ESIO_LIB) -lesiof -lm
ESIO_INCLUDE = -I$(ESIO_LIB)
HDF_LIB=-L$(HDF_DIR)/lib -lhdf5_hl -lhdf5 -lm -lz
HDF_INC=-I$(HDF_DIR)/include
endif
########################################################################


LIBS=$(GSL_LIBRARY) $(FFTW3_LIB)  $(ESIO_LIBRARY) $(HDF_LIB) $(OTHERLIBS)
INCLUDE= $(FFTW3_INC) $(GSL_INCLUDE) $(HDF_INC) $(ESIO_INCLUDE) -I$(ESIO_DIR)/lib $(OTHERINC)

## build the executable ##
$(EXEC): $(OBJECTS)
	$(COMPILE) $(INCLUDE) $(OBJECTS) -o  $@ $(LDFLAGS) $(LIBS) 

# build cpp object files
%.o: %.cpp
	@echo building $< 
	$(CC) -c $(CFLAGS) $(INCLUDE) $< -o $@

# build f90 object files
%.o: %.f90
	@echo building $< 
	$(FC) -c $(CFLAGS) $(INCLUDE) $< -o $@

# build f object files
%.o: %.f
	@echo building $< 
	$(FC) -c $(CFLAGS) $(INCLUDE) $< -o $@


# build c object files
%.o: %.c
	@echo building $<
	$(CC) -c $(CFLAGS) $(INCLUDE) $< -o $@

# clean directive 'make clean'
clean:
	- /bin/rm -f $(EXEC) *.o *.mod *~ \#* 
	@echo 'files cleaned'

