.SUFFIXES: .F90 .o
  
all: dummy MPAS_SMOKE

dummy:
	echo "****** compiling MPAS_SMOKE ******"

OBJS = \
	mpas_smoke_config.o  \
	dep_dry_simple_mod.o \
	dep_dry_mod_emerson.o \
	dep_data_mod.o \
	module_wetdep_ls.o \
	module_zero_plumegen_coms.o \
	dust_data_mod.o \
	dust_fengsha_mod.o \
	ssalt_mod.o \
	seas_ngac_mod.o \
	seas_data_mod.o \
	module_smoke_plumerise.o \
	module_add_emiss_burn.o \
	module_plumerise.o \
	mpas_smoke_wrapper.o \
	mpas_smoke_init.o
   

MPAS_SMOKE: $(OBJS)
	ar -ru ./../../libphys.a $(OBJS)

# DEPENDENCIES:
#
mpas_smoke_init.o:

dep_dry_mod_emerson.o: \
        dep_data_mod.o \
	mpas_smoke_init.o \
	mpas_smoke_config.o

dust_fengsha_mod.o: \
	mpas_smoke_config.o \
	dust_data_mod.o \
	mpas_smoke_init.o

dust_data_mod.o: \
	mpas_smoke_init.o

module_wetdep_ls.o: \
	mpas_smoke_config.o \
	mpas_smoke_init.o

module_add_emiss_burn.o: \
	mpas_smoke_config.o

module_smoke_plumerise.o:  \
	module_zero_plumegen_coms.o

module_plumerise.o: \
	module_smoke_plumerise.o \
	module_zero_plumegen_coms.o \
	mpas_smoke_config.o

mpas_smoke_wrapper.o:      \
        mpas_smoke_config.o \
        module_add_emiss_burn.o \
        module_plumerise.o \
	module_wetdep_ls.o \
        dep_dry_simple_mod.o \
        dep_dry_mod_emerson.o

module_zero_plumegen_coms.o: \ \
        mpas_smoke_config.o

ssalt_mod.o: \
	seas_ngac_mod.o \
	seas_data_mod.o \
	mpas_smoke_init.o

seas_data_mod.o:

seas_ngac_mod.o:

#module_smoke_plumerise.o: module_smoke_plumerise.F90
#	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS) -c module_smoke_plumerise.F90 $(CPPINCLUDES) $(FCINCLUDES) -I.. -I../../../framework -I../../../external/esmf_time_f90

#module_zero_plumegen_coms.o: module_zero_plumegen_coms.F90
#	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS) -c module_zero_plumegen_coms.F90 $(CPPINCLUDES) $(FCINCLUDES) -I.. -I../../../framework -I../../../external/esmf_time_f90

clean:
	$(RM) *.f90 *.o *.mod
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.F90.o:
ifeq "$(GEN_F90)" "true"
        $(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) $< > $*.f90
        $(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I.. -I../../../../framework -I../../../../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS) -c $*.F90 $(CPPINCLUDES) $(FCINCLUDES) -I.. -I../../../../framework -I../../../../external/esmf_time_f90
endif
