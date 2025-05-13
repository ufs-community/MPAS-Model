.SUFFIXES: .F90 .o

all: dummy lsm_ruc

dummy:
	echo "****** compiling RUC land and ice models ******"

OBJS = \
	module_ruc_land.o	\
	module_ruc_ice.o		\
	module_sf_sfcdiags_ruclsm.o

lsm_ruc: $(OBJS)
	ar -ru ./../../libphys.a $(OBJS)

# DEPENDENCIES:
module_ruc_land.o: \
        ../../mpas_atmphys_constants.o
module_ruc_ice.o: \
        ../../mpas_atmphys_constants.o
#

clean:
	$(RM) *.F90 *.F *.o *.mod
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.F90.o:
ifeq "$(GEN_F90)" "true"
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../.. -I../../physics_wrf -I../../physics_mmm -I../../../../framework -I../../../../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS) -c $*.F90 $(CPPINCLUDES) $(FCINCLUDES) -I../.. -I../../physics_wrf -I../../physics_mmm -I../../../../framework -I../../../../external/esmf_time_f90
endif
