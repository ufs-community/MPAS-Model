# Every physics directory that holds Fortran module files. Some compilers (nvfortran) also need
# the module files of modules that the physics modules used here use in turn.
PHYSICS_MODDIRS = $(addprefix -I../physics/, physics_mmm physics_wrf physics_noahmp physics_noahmp/utility physics_noahmp/drivers/mpas physics_noahmp/src physics_noaa/TEMPO physics_noaa/UGWP physics_noaa/MYNN-EDMF/src physics_noaa/MYNN-SFC/src physics_noaa/RUCLSM physics_noaa/GFL physics_noaa/SMOKE)

.SUFFIXES: .F .o

OBJS = mpas_atm_time_integration.o \
       mpas_atm_boundaries.o \
       mpas_atm_dissipation_models.o

all: $(OBJS)

mpas_atm_time_integration.o: mpas_atm_boundaries.o mpas_atm_iau.o mpas_atm_dissipation_models.o

mpas_atm_boundaries.o:


clean:
	$(RM) *.o *.mod *.f90
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(PHYSICS) $(CPPINCLUDES) -I../../framework $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I.. -I../../framework -I../../operators -I../physics -I../physics/physics_wrf -I../physics/physics_mmm $(MPAS_ESMF_INC) $(PHYSICS_MODDIRS)
else
	$(FC) $(CPPFLAGS) $(PHYSICS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I.. -I../../framework -I../../operators -I../physics -I../physics/physics_wrf -I../physics/physics_mmm $(MPAS_ESMF_INC) -I../stochastic_physics $(PHYSICS_MODDIRS)
endif
