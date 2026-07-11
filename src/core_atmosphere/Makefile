.SUFFIXES: .F .o

#
# To build a dycore-only MPAS-Atmosphere model, comment-out or delete
# the definition of PHYSICS and CHEMISTRY, below
#
# If MPAS_CAM_DYCORE is found in CPPFLAGS, PHYSICS and CHEMISTRY will become undefined automatically
#
ifeq ($(findstring MPAS_CAM_DYCORE,$(CPPFLAGS)),)
    PHYSICS = -DDO_PHYSICS
    CHEMISTRY = -DDO_CHEMISTRY
endif

ifdef PHYSICS
    PHYSCORE = physcore
    PHYS_OBJS = libphys/*.o
endif

ifdef CHEMISTRY
    CHEMCORE = chemcore
    CHEM_OBJS = libchem/*.o
endif

OBJS = mpas_atm_core.o \
       mpas_atm_core_interface.o \
       mpas_atm_dimensions.o \
       mpas_atm_threading.o \
       mpas_atm_halos.o

all: $(PHYSCORE) $(CHEMCORE) dycore diagcore atmcore utilities

core_reg:
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $(PHYSICS) $(CHEMISTRY) Registry.xml > Registry_processed.xml

core_input_gen:
	if [ ! -e default_inputs ]; then  mkdir default_inputs; fi
	rm -f default_inputs/streams.atmosphere default_inputs/stream_list.atmosphere.*
	( cd default_inputs; $(NL_GEN) ../Registry_processed.xml namelist.atmosphere in_defaults=true )
	( cd default_inputs; $(ST_GEN) ../Registry_processed.xml streams.atmosphere stream_list.atmosphere. listed in_defaults=true)

gen_includes: core_reg
	(if [ ! -d inc ]; then mkdir -p inc; fi) # To generate *.inc files
	(cd inc; $(REG_PARSE) ../Registry_processed.xml $(CPPFLAGS) )

post_build:
	if [ ! -e $(ROOT_DIR)/default_inputs ]; then mkdir $(ROOT_DIR)/default_inputs; fi
	rm -f $(ROOT_DIR)/default_inputs/streams.$(CORE) $(ROOT_DIR)/default_inputs/stream_list.$(CORE).*
	cp default_inputs/* $(ROOT_DIR)/default_inputs/.
	( cd $(ROOT_DIR)/default_inputs; for FILE in `ls -1`; do \
	    if [ -e ../$$FILE ]; then \
	        N=1; \
	        while [ -e "../$$FILE.old$$(printf '%03d' $$N)" ]; do N=$$((N+1)); done; \
	        mv "../$$FILE" "../$$FILE.old$$(printf '%03d' $$N)"; \
	    fi; \
	    cp $$FILE ../.; \
	done )

physcore: mpas_atm_dimensions.o
	( cd physics; $(MAKE) all )
	( cd stochastic_physics; $(MAKE) -f Makefile all )
	( mkdir libphys; cd libphys; ar -x ../physics/libphys.a; ar -x ../stochastic_physics/libstochphys.a ) 
	( cd ../..; ln -sf ./src/core_atmosphere/physics/physics_wrf/files/*TBL .)
	( cd ../..; ln -sf ./src/core_atmosphere/physics/physics_wrf/files/*DATA* .)
	( cd ../..; ln -sf ./src/core_atmosphere/physics/physics_noahmp/parameters/*TBL .)

chemcore:
	( cd chemistry; $(MAKE) all CHEMISTRY="$(CHEMISTRY)" )
	( mkdir libchem; cd libchem; ar -x ../chemistry/libchem.a )

dycore: mpas_atm_dimensions.o $(PHYSCORE) $(CHEMCORE)
	( cd dynamics; $(MAKE) all PHYSICS="$(PHYSICS)" )

diagcore: $(PHYSCORE) dycore
	( cd diagnostics; $(MAKE) all PHYSICS="$(PHYSICS)" )

utilities: $(PHYSCORE)
	( cd utils; $(MAKE) all PHYSICS="$(PHYSICS)" )

atmcore: $(PHYSCORE) dycore diagcore $(OBJS)
	ar -ru libdycore.a $(OBJS) dynamics/*.o $(PHYS_OBJS) $(CHEM_OBJS) diagnostics/*.o

mpas_atm_core_interface.o: mpas_atm_core.o

mpas_atm_core.o: dycore diagcore mpas_atm_threading.o mpas_atm_halos.o

mpas_atm_dimensions.o:

clean:
	( cd physics; $(MAKE) clean )
	( cd chemistry; $(MAKE) clean )
	( cd dynamics; $(MAKE) clean )
	@# The stochastic_physics directory contains both a makefile (lower case
	@# "m") and a Makefile (upper case "M").  The former is used for FV3 and
	@# the latter for MPAS.  The default for "make" is to use "makefile", but
	@# here, we want to use the one for MPAS (with an upper case "M").  Indicate
	@# that using the -f flag.
	( cd stochastic_physics; $(MAKE) -f Makefile clean )
	( cd diagnostics; $(MAKE) clean )
	( cd utils; $(MAKE) clean )
	( cd ../..; rm -f *TBL )
	( cd ../..; rm -f *DATA* )
	$(RM) -r libphys
	$(RM) *.o *.mod *.f90 libdycore.a
	$(RM) Registry_processed.xml
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i
	$(RM) -r default_inputs

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(PHYSICS) $(CHEMISTRY) $(CPPINCLUDES) -I./inc -I../framework $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework -I../operators -I./physics -I./dynamics -I./diagnostics -I./physics/physics_wrf -I./physics/physics_mmm -I./physics/physics_noaa/UGWP $(MPAS_ESMF_INC) -I./chemistry -I./stochastic_physics
else
	$(FC) $(CPPFLAGS) $(PHYSICS) $(CHEMISTRY) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I./inc -I../framework -I../operators -I./physics -I./dynamics -I./diagnostics -I./physics/physics_wrf -I./physics/physics_mmm -I./physics/physics_noaa/UGWP $(MPAS_ESMF_INC) -I./chemistry -I./stochastic_physics
endif
