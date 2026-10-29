.SUFFIXES: .F .f .f90 .o

OBJS = \
	gfs_tke_edmf_compat.o \
	tridi_mod.o           \
	mfpbltq_mod.o         \
	mfscuq_mod.o          \
	canopy_utils_mod.o    \
	satmedmfvdifq.o       \
	progsigma_calc.o     \
	progomega_calc.o     \
	samfaerosols.o        \
	samfdeepcnv.o         \
	samfshalcnv.o

all: $(OBJS)
	ar -ru libphys_gfs.a $(OBJS)

physics_gfs_lib: $(OBJS)
	ar -ru ../libphys.a $(OBJS)

gfs_tke_edmf_compat.o: \
	../mpas_atmphys_constants.o

tridi_mod.o: \
	gfs_tke_edmf_compat.o

mfpbltq_mod.o: \
	gfs_tke_edmf_compat.o

mfscuq_mod.o: \
	gfs_tke_edmf_compat.o

canopy_utils_mod.o: \
	gfs_tke_edmf_compat.o

satmedmfvdifq.o: \
	gfs_tke_edmf_compat.o \
	tridi_mod.o \
	mfpbltq_mod.o \
	mfscuq_mod.o \
	canopy_utils_mod.o

progsigma_calc.o: \
	gfs_tke_edmf_compat.o

progomega_calc.o: \
	gfs_tke_edmf_compat.o

samfaerosols.o: \
	gfs_tke_edmf_compat.o

samfdeepcnv.o: \
	gfs_tke_edmf_compat.o \
	progsigma_calc.o \
	progomega_calc.o \
	samfaerosols.o

samfshalcnv.o: \
	gfs_tke_edmf_compat.o \
	progsigma_calc.o \
	progomega_calc.o \
	samfaerosols.o

# Select source-form flags for supported compiler families.
# Keep these flags after FFLAGS so they override the global source form.
# BUILD_TARGET identifies the toolchain even when FC is an MPI wrapper.
GFS_COMPILER = $(notdir $(firstword $(FC)))
ifneq (,$(filter nvhpc pgi,$(BUILD_TARGET)))
  GFS_FIXED_FLAGS = -Mfixed -Mextend
  GFS_FREE_FLAGS  = -Mfree
else ifneq (,$(findstring nvhpc,$(BUILD_TARGET)))
  GFS_FIXED_FLAGS = -Mfixed -Mextend
  GFS_FREE_FLAGS  = -Mfree
else ifneq (,$(findstring pgi,$(BUILD_TARGET)))
  GFS_FIXED_FLAGS = -Mfixed -Mextend
  GFS_FREE_FLAGS  = -Mfree
else ifneq (,$(filter nvfortran pgf90 mpipgf90,$(GFS_COMPILER)))
  GFS_FIXED_FLAGS = -Mfixed -Mextend
  GFS_FREE_FLAGS  = -Mfree
else ifneq (,$(findstring intel,$(BUILD_TARGET)))
  GFS_FIXED_FLAGS = -fixed -extend-source 132
  GFS_FREE_FLAGS  = -free
else ifneq (,$(filter ifort ifx mpiifort mpiifx,$(GFS_COMPILER)))
  GFS_FIXED_FLAGS = -fixed -extend-source 132
  GFS_FREE_FLAGS  = -free
else ifneq (,$(findstring gnu,$(BUILD_TARGET)))
  GFS_FIXED_FLAGS = -ffixed-form -ffixed-line-length-none
  GFS_FREE_FLAGS  = -ffree-form -ffree-line-length-none
else ifneq (,$(filter gfortran mpif90-gfortran,$(GFS_COMPILER)))
  GFS_FIXED_FLAGS = -ffixed-form -ffixed-line-length-none
  GFS_FREE_FLAGS  = -ffree-form -ffree-line-length-none
else
  $(error Unsupported GFS Fortran compiler: BUILD_TARGET=$(BUILD_TARGET) FC=$(FC))
endif

GFS_FIXED_OBJECTS = gfs_tke_edmf_compat.o mfscuq_mod.o \
                    canopy_utils_mod.o samfaerosols.o
$(GFS_FIXED_OBJECTS): GFS_SOURCE_FORM_FLAGS = $(GFS_FIXED_FLAGS)

GFS_FREE_F_OBJECTS = samfdeepcnv.o samfshalcnv.o
$(GFS_FREE_F_OBJECTS): GFS_SOURCE_FORM_FLAGS = $(GFS_FREE_FLAGS)

clean:
	$(RM) *.o *.mod *.i libphys_gfs.a
	$(RM) -f gfs_tke_edmf_compat.f90 tridi_mod.f90 mfpbltq_mod.f90 mfscuq_mod.f90 canopy_utils_mod.f90 satmedmfvdifq.f90 samfaerosols.f90

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(COREDEF) $(HYDROSTATIC) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) $(GFS_SOURCE_FORM_FLAGS) -c $*.f90 $(FCINCLUDES) -I. -I.. -I../physics_mmm -I../physics_wrf -I../../../framework -I../../../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(HYDROSTATIC) $(FFLAGS) $(GFS_SOURCE_FORM_FLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I. -I.. -I../physics_mmm -I../physics_wrf -I../../../framework -I../../../external/esmf_time_f90
endif

.f.o:
	$(RM) $@ $*.mod
	$(FC) $(FFLAGS) $(GFS_SOURCE_FORM_FLAGS) -c $< $(FCINCLUDES) -I. -I.. -I../physics_mmm -I../physics_wrf -I../../../framework -I../../../external/esmf_time_f90

.f90.o:
	$(RM) $@ $*.mod
	$(FC) $(FFLAGS) $(GFS_SOURCE_FORM_FLAGS) -c $< $(FCINCLUDES) -I. -I.. -I../physics_mmm -I../physics_wrf -I../../../framework -I../../../external/esmf_time_f90
