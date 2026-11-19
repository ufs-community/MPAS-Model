# Every physics directory that holds Fortran module files. Some compilers (nvfortran) also need
# the module files of modules that the physics modules used here use in turn.
PHYSICS_MODDIRS = $(addprefix -I../physics/, physics_mmm physics_wrf physics_noahmp physics_noahmp/utility physics_noahmp/drivers/mpas physics_noahmp/src physics_noaa/TEMPO physics_noaa/UGWP physics_noaa/MYNN-EDMF/src physics_noaa/MYNN-SFC/src physics_noaa/RUCLSM physics_noaa/GFL physics_noaa/SMOKE)

.SUFFIXES: .F .o

#
# Add new diagnostic modules to DIAGNOSTIC_MODULES
#
DIAGNOSTIC_MODULES = \
	mpas_atm_diagnostic_template.o \
	mpas_isobaric_diagnostics.o \
	mpas_cloud_diagnostics.o \
	mpas_convective_diagnostics.o \
	mpas_pv_diagnostics.o \
	mpas_soundings.o \

mpas_isobaric_diagnostics.o: mpas_atm_diagnostics_utils.o

mpas_cloud_diagnostics.o: mpas_atm_diagnostics_utils.o

mpas_convective_diagnostics.o: mpas_atm_diagnostics_utils.o

mpas_pv_diagnostics.o: mpas_atm_diagnostics_utils.o

mpas_soundings.o:


################### Generally no need to modify below here ###################


OBJS = mpas_atm_diagnostics_manager.o mpas_atm_diagnostics_utils.o

all: $(DIAGNOSTIC_MODULES) $(OBJS)

mpas_atm_diagnostics_manager.o: mpas_atm_diagnostics_utils.o $(DIAGNOSTIC_MODULES)


clean:
	$(RM) *.o *.mod *.f90
	@# Some Intel compilers generate *.i files; clean them up, too
	$(RM) *.i

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(PHYSICS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../../framework -I../../operators -I../dynamics -I../physics -I../physics/physics_wrf $(MPAS_ESMF_INC) $(PHYSICS_MODDIRS)
else
	$(FC) $(CPPFLAGS) $(PHYSICS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../../framework -I../../operators -I../dynamics -I../physics -I../physics/physics_wrf $(MPAS_ESMF_INC) $(PHYSICS_MODDIRS)
endif
