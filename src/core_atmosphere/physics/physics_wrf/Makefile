.SUFFIXES: .F .o .F90

.PHONY: physics_wrf physics_wrf_lib

all: dummy physics_wrf

dummy:
	echo "****** compiling physics_wrf ******"

OBJS = \
	libmassv.o                     \
	module_bep_bem_helper.o        \
	module_bl_gwdo.o               \
	module_bl_ugwp_gwdo.o          \
	module_bl_mynn.o               \
	module_bl_ysu.o                \
	module_bl_myjpbl.o             \
	module_cam_error_function.o    \
	module_cam_shr_kind_mod.o      \
	module_cam_support.o           \
	module_cu_gf.mpas.o            \
	module_cu_tiedtke.o            \
	module_cu_ntiedtke.o           \
	module_cu_kfeta.o              \
	module_mp_kessler.o            \
	module_mp_nssl_2mom.o          \
	module_mp_radar.o              \
	module_mp_thompson.o           \
	module_mp_thompson_aerosols.o  \
	module_mp_thompson_cldfra3.o   \
	module_mp_wsm6.o               \
	module_ra_cam.o                \
	module_ra_cam_support.o        \
	module_ra_rrtmg_lw.o           \
	module_ra_rrtmg_sw.o           \
	module_ra_rrtmg_sw_aerosols.o  \
	module_ra_rrtmg_vinterp.o      \
	module_sf_bem.o                \
	module_sf_bep.o                \
	module_sf_bep_bem.o            \
	module_sf_sfcdiags.o           \
	module_sf_myjsfc.o             \
	module_sf_mynn.o               \
	module_sf_noahdrv.o            \
	module_sf_noahlsm.o            \
        module_sf_noahlsm_glacial_only.o \
        module_sf_noah_seaice.o        \
        module_sf_noah_seaice_drv.o    \
	module_sf_oml.o                \
	module_sf_sfclay.o             \
	module_sf_sfclayrev.o          \
	module_sf_urban.o              \
	bl_mynn_post.o                 \
	bl_mynn_pre.o                  \
	module_mynnsfclay.o	       \
        module_mynnsfclay_pre.o	       \
	module_mynnsfclay_driver.o     \
	cu_ntiedtke_post.o             \
	cu_ntiedtke_pre.o	       \
	sf_mynn_pre.o                  \
	sf_sfclayrev_pre.o


physics_wrf: $(OBJS)

physics_wrf_lib:
	ar -ru ./../libphys.a $(OBJS)

# DEPENDENCIES:
module_bl_mynn.o: \
        bl_mynn_post.o \
        bl_mynn_pre.o

module_mynnsfclay_driver.o: \
	module_mynnsfclay.o \
	module_mynnsfclay_pre.o

module_cam_support.o: \
	module_cam_shr_kind_mod.o

module_cu_ntiedtke.o: \
	cu_ntiedtke_post.o \
	cu_ntiedtke_pre.o

module_ra_cam.o: \
	module_cam_support.o \
	module_ra_cam_support.o

module_ra_rrtmg_lw.o: \
	module_ra_rrtmg_vinterp.o

module_ra_rrtmg_sw.o: \
	module_ra_rrtmg_lw.o \
	module_ra_rrtmg_vinterp.o

module_sf_bep.o: \
	module_bep_bem_helper.o \
	module_sf_urban.o

module_sf_bep_bem.o: \
	module_bep_bem_helper.o \
	module_sf_bem.o \
	module_sf_urban.o

module_sf_mynn.o: \
	sf_mynn_pre.o

module_sf_sfclayrev.o: \
	sf_sfclayrev_pre.o

module_sf_noahdrv.o: \
	module_sf_bem.o \
	module_sf_bep.o \
	module_sf_bep_bem.o \
	module_sf_noahlsm.o \
	module_sf_noahlsm_glacial_only.o \
	module_sf_urban.o

module_sf_noahlsm_glacial_only.o: \
	module_sf_noahlsm.o

module_sf_noah_seaice_drv.o: \
	module_sf_noah_seaice.o

module_sf_noah_seaice.o: \
	module_sf_noahlsm.o

clean:
	$(RM) *.f90 *.o *.mod
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.F.o:
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(COREDEF) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I.. -I../physics_mmm -I../physics_noaa/UGWP -I../../../framework -I../../../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I.. -I../physics_mmm -I../physics_noaa/UGWP -I../../../framework -I../../../external/esmf_time_f90
endif

.F90.o:
ifeq "$(GEN_F90)" "true"
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I.. -I../physics_mmm -I../physics_noaa/UGWP -I../../../framework -I../../../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS) -c $*.F90 $(CPPINCLUDES) $(FCINCLUDES) -I.. -I../physics_mmm -I../physics_noaa/UGWP -I../../../framework -I../../../external/esmf_time_f90
endif
