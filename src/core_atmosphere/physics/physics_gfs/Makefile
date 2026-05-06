.SUFFIXES: .F .o

OBJS = \
	gfs_tke_edmf_compat.o \
	tridi_mod.o           \
	mfpbltq_mod.o         \
	mfscuq_mod.o          \
	canopy_utils_mod.o    \
	satmedmfvdifq.o

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

clean:
	$(RM) *.o *.mod *.f90 *.i libphys_gfs.a

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(COREDEF) $(HYDROSTATIC) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I. -I.. -I../physics_mmm -I../physics_wrf -I../../../framework -I../../../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(COREDEF) $(HYDROSATIC) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I. -I.. -I../physics_mmm -I../physics_wrf -I../../../framework -I../../../external/esmf_time_f90
endif
