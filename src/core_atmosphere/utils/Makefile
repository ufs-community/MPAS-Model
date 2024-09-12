.SUFFIXES: .F .o

ifdef PHYSICS
    UTILS = build_tables
endif

all: $(UTILS)

build_tables: build_tables.o atmphys_build_tables_thompson.o
	$(LINKER) $(LDFLAGS) -o build_tables build_tables.o atmphys_build_tables_thompson.o -L../../framework -L../physics -lphys -lframework $(LIBS) -L../../external/esmf_time_f90 -lesmf_time -L/apps/netcdf/4.2.1.1-intel/lib -lnetcdff -I/apps/netcdf/4.2.1.1-intel/include
	mv build_tables ../../..


build_tables.o: \
	atmphys_build_tables_thompson.o

atmphys_build_tables_thompson.o: \
	../physics/TEMPO/module_mp_thompson.o \
	../physics/TEMPO/module_mp_thompson_utils.o \
	../physics/TEMPO/module_mp_thompson_params.o

clean:
	$(RM) ../../../build_tables
	$(RM) *.o *.mod *.f90
	@# Some Intel compilers generate *.i files; clean them up, too
	$(RM) *.i

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(PHYSICS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../../framework -I../../operators -I../physics -I../physics/physics_mmm -I../physics/physics_wrf -I../physics/TEMPO -I../../external/esmf_time_f90 -L/apps/netcdf/4.2.1.1-intel/lib -lnetcdff -I/apps/netcdf/4.2.1.1-intel/include
else
	$(FC) $(CPPFLAGS) $(PHYSICS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../../framework -I../../operators -I../physics -I../physics/physics_mmm -I../physics/physics_wrf -I../physics/TEMPO -I../../external/esmf_time_f90 -L/apps/netcdf/4.2.1.1-intel/lib -lnetcdff -I/apps/netcdf/4.2.1.1-intel/include
endif
