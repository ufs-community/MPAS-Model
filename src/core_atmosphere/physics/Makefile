.SUFFIXES: .F .o .F90

ifeq ($(CORE),atmosphere)
COREDEF = -Dmpas
endif

all: lookup_tables core_physics_init core_physics_mmm core_microphysics core_UGWP core_physics_wrf core_physics

dummy:
	echo "****** compiling physics ******"

OBJS_init = \
	ccpp_kinds.o             \
	mpas_atmphys_constants.o \
	mpas_atmphys_date_time.o \
	mpas_atmphys_functions.o \
	mpas_atmphys_utilities.o

OBJS = \
	mpas_atmphys_camrad_init.o         \
	mpas_atmphys_control.o             \
	mpas_atmphys_driver.o              \
	mpas_atmphys_driver_cloudiness.o   \
	mpas_atmphys_driver_convection.o   \
	mpas_atmphys_driver_gwdo.o         \
	mpas_atmphys_driver_lsm.o          \
	mpas_atmphys_driver_microphysics.o \
	mpas_atmphys_driver_oml.o          \
	mpas_atmphys_driver_pbl.o          \
	mpas_atmphys_driver_radiation_lw.o \
	mpas_atmphys_driver_radiation_sw.o \
	mpas_atmphys_driver_seaice.o       \
	mpas_atmphys_driver_sfclayer.o     \
	mpas_atmphys_finalize.o            \
	mpas_atmphys_init.o                \
	mpas_atmphys_init_microphysics.o   \
	mpas_atmphys_interface.o           \
	mpas_atmphys_landuse.o             \
	mpas_atmphys_lsm_noahinit.o        \
	mpas_atmphys_lsm_shared.o          \
	mpas_atmphys_lsm_rucinit.o         \
	mpas_atmphys_manager.o             \
	mpas_atmphys_o3climatology.o       \
	mpas_atmphys_packages.o            \
	mpas_atmphys_rrtmg_lwinit.o        \
	mpas_atmphys_rrtmg_swinit.o        \
	mpas_atmphys_todynamics.o          \
	mpas_atmphys_update_surface.o      \
	mpas_atmphys_update.o              \
	mpas_atmphys_vars.o

lookup_tables:
	./checkout_data_files.sh

core_physics_mmm: core_physics_init
	(cd physics_mmm; $(MAKE) all)

core_microphysics: core_physics_init core_physics_mmm
	(if [ -f physics_wrf/module_mp_thompson.F ]; then rm physics_wrf/module_mp_thompson.F; else echo "No old module_mp_thompson.F file found"; fi)
	(cd TEMPO; cp ./drivers/mpas/module_mp_thompson.F90 .; $(MAKE) all COREDEF="$(COREDEF)")

core_UGWP: core_physics_init
	(cd UGWP; $(MAKE) all)

core_physics_wrf: core_physics_init core_physics_mmm core_UGWP
	(cd physics_wrf; $(MAKE) all COREDEF="$(COREDEF)")

core_physics_init: $(OBJS_init)
	ar -ru libphys.a $(OBJS_init)

core_physics: core_microphysics core_physics_wrf
	($(MAKE) phys_interface COREDEF="$(COREDEF)")
	ar -ru libphys.a $(OBJS)

phys_interface: $(OBJS)

# DEPENDENCIES:
mpas_atmphys_camrad_init.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_utilities.o

mpas_atmphys_control.o: \
	mpas_atmphys_utilities.o \
	mpas_atmphys_vars.o

mpas_atmphys_driver.o: \
	mpas_atmphys_driver_cloudiness.o \
	mpas_atmphys_driver_convection.o \
	mpas_atmphys_driver_gwdo.o \
	mpas_atmphys_driver_lsm.o \
	mpas_atmphys_driver_pbl.o \
	mpas_atmphys_driver_radiation_lw.o \
	mpas_atmphys_driver_radiation_sw.o \
	mpas_atmphys_driver_seaice.o \
	mpas_atmphys_driver_sfclayer.o \
	mpas_atmphys_driver_oml.o \
	mpas_atmphys_constants.o \
	mpas_atmphys_interface.o \
	mpas_atmphys_update.o \
	mpas_atmphys_vars.o

mpas_atmphys_driver_cloudiness.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_vars.o

mpas_atmphys_driver_convection.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_utilities.o \
	mpas_atmphys_vars.o

mpas_atmphys_driver_gwdo.o: \
	mpas_atmphys_vars.o

mpas_atmphys_driver_lsm.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_landuse.o \
	mpas_atmphys_lsm_noahinit.o \
	mpas_atmphys_lsm_rucinit.o \
	mpas_atmphys_lsm_shared.o \
	mpas_atmphys_vars.o 

mpas_atmphys_driver_microphysics.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_init_microphysics.o \
	mpas_atmphys_interface.o \
	mpas_atmphys_vars.o

mpas_atmphys_driver_oml.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_landuse.o \
	mpas_atmphys_vars.o

mpas_atmphys_driver_pbl.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_vars.o

mpas_atmphys_driver_radiation_lw.o: \
	mpas_atmphys_camrad_init.o \
	mpas_atmphys_constants.o \
	mpas_atmphys_driver_radiation_sw.o \
	mpas_atmphys_manager.o \
	mpas_atmphys_rrtmg_lwinit.o \
	mpas_atmphys_vars.o

mpas_atmphys_driver_radiation_sw.o: \
	mpas_atmphys_camrad_init.o \
	mpas_atmphys_constants.o \
	mpas_atmphys_manager.o \
	mpas_atmphys_rrtmg_swinit.o \
	mpas_atmphys_vars.o

mpas_atmphys_driver_sfclayer.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_vars.o

mpas_atmphys_init.o: \
	mpas_atmphys_driver_convection.o \
	mpas_atmphys_driver_lsm.o \
	mpas_atmphys_driver_microphysics.o \
	mpas_atmphys_driver_pbl.o \
	mpas_atmphys_driver_radiation_lw.o \
	mpas_atmphys_driver_radiation_sw.o \
	mpas_atmphys_driver_sfclayer.o \
	mpas_atmphys_landuse.o \
	mpas_atmphys_o3climatology.o \
	mpas_atmphys_vars.o

mpas_atmphys_interface.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_vars.o

mpas_atmphys_landuse.o: \
	mpas_atmphys_utilities.o \
	mpas_atmphys_vars.o

mpas_atmphys_lsm_noahinit.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_utilities.o

mpas_atmphys_lsm_rucinit.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_utilities.o

mpas_atmphys_manager.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_o3climatology.o \
	mpas_atmphys_update_surface.o \
	mpas_atmphys_utilities.o \
	mpas_atmphys_vars.o

mpas_atmphys_o3climatology.o: \
	mpas_atmphys_date_time.o \
	mpas_atmphys_constants.o \
	mpas_atmphys_utilities.o

mpas_atmphys_rrtmg_lwinit.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_utilities.o

mpas_atmphys_rrtmg_swinit.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_utilities.o

mpas_atmphys_driver_seaice.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_lsm_shared.o \
	mpas_atmphys_vars.o

mpas_atmphys_todynamics.o: \
	mpas_atmphys_constants.o \
	mpas_atmphys_vars.o

mpas_atmphys_update_surface.o: \
	mpas_atmphys_date_time.o \
	mpas_atmphys_constants.o \
	mpas_atmphys_vars.o

mpas_atmphys_update.o: \
	mpas_atmphys_driver_convection.o \
	mpas_atmphys_vars.o

clean:
	$(RM) *.o *.mod *.f90 libphys.a
	( cd physics_wrf; $(MAKE) clean )
	( cd physics_mmm; $(MAKE) clean )
	( cd TEMPO; $(MAKE) clean )
	( cd UGWP; $(MAKE) clean )
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(COREDEF) $(HYDROSTATIC) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I./physics_mmm -I./physics_wrf -I./TEMPO -I./UGWP -I.. -I../../framework -I../../external/esmf_time_f90 -L/apps/netcdf/4.2.1.1-intel/lib -lnetcdff -I/apps/netcdf/4.2.1.1-intel/include -L/apps/netcdf/4.2.1.1-intel/lib -lnetcdff -I/apps/netcdf/4.2.1.1-intel/include
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(HYDROSATIC) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I./physics_mmm -I./physics_wrf -I./TEMPO -I./UGWP -I.. -I../../framework -I../../external/esmf_time_f90 -L/apps/netcdf/4.2.1.1-intel/lib -lnetcdff -I/apps/netcdf/4.2.1.1-intel/include -L/apps/netcdf/4.2.1.1-intel/lib -lnetcdff -I/apps/netcdf/4.2.1.1-intel/include
endif

.F90.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(COREDEF) $(HYDROSTATIC) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I./physics_mmm -I./physics_wrf -I./TEMPO -I./UGWP -I.. -I../../framework -I../../external/esmf_time_f90 -L/apps/netcdf/4.2.1.1-intel/lib -lnetcdff -I/apps/netcdf/4.2.1.1-intel/include -L/apps/netcdf/4.2.1.1-intel/lib -lnetcdff -I/apps/netcdf/4.2.1.1-intel/include
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(HYDROSATIC) $(FFLAGS) -c $*.F90 $(CPPINCLUDES) $(FCINCLUDES) -I./physics_mmm -I./physics_wrf -I./TEMPO -I./UGWP -I.. -I../../framework -I../../external/esmf_time_f90 -L/apps/netcdf/4.2.1.1-intel/lib -lnetcdff -I/apps/netcdf/4.2.1.1-intel/include -L/apps/netcdf/4.2.1.1-intel/lib -lnetcdff -I/apps/netcdf/4.2.1.1-intel/include
endif
