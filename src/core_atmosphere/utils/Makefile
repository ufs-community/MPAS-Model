# Every physics directory that holds Fortran module files. Some compilers (nvfortran) also need
# the module files of modules that the physics modules used here use in turn.
PHYSICS_MODDIRS = $(addprefix -I../physics/, physics_mmm physics_wrf physics_noahmp physics_noahmp/utility physics_noahmp/drivers/mpas physics_noahmp/src physics_noaa/TEMPO physics_noaa/UGWP physics_noaa/MYNN-EDMF/src physics_noaa/MYNN-SFC/src physics_noaa/RUCLSM physics_noaa/GFL physics_noaa/SMOKE)

.SUFFIXES: .F .o

ifdef PHYSICS
    UTILS = build_tables
endif

all: $(UTILS)

build_tables: build_tables.o atmphys_build_tables_thompson.o
	$(LINKER) $(LDFLAGS) -o build_tables build_tables.o atmphys_build_tables_thompson.o -L../../framework -L../physics -lphys -lframework $(LIBS) $(MPAS_ESMF_LIB)
	mv build_tables ../../..

build_tables.o: \
	atmphys_build_tables_thompson.o

clean:
	$(RM) ../../../build_tables
	$(RM) *.o *.mod *.f90
	@# Some Intel compilers generate *.i files; clean them up, too
	$(RM) *.i

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(PHYSICS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../../framework -I../../operators -I../physics -I../physics/physics_mmm -I../physics/physics_wrf $(MPAS_ESMF_INC) $(PHYSICS_MODDIRS)
else
	$(FC) $(CPPFLAGS) $(PHYSICS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../../framework -I../../operators -I../physics -I../physics/physics_mmm -I../physics/physics_wrf $(MPAS_ESMF_INC) $(PHYSICS_MODDIRS)
endif
