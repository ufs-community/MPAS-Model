.SUFFIXES: .F .o

#
# To build a dycore-only MPAS-Atmosphere model, comment-out or delete
# the definition of PHYSICS, below
#
# If MPAS_CAM_DYCORE is found in CPPFLAGS, PHYSICS will become undefined automatically
#
ifeq ($(findstring MPAS_CAM_DYCORE,$(CPPFLAGS)),)
    PHYSICS = -DDO_PHYSICS
endif

ifdef PHYSICS
    PHYSCORE = physcore
    PHYS_OBJS = libphys/*.o
endif

OBJS = mpas_atm_core.o \
       mpas_atm_core_interface.o \
       mpas_atm_dimensions.o \
       mpas_atm_threading.o \
       mpas_atm_halos.o

all: $(PHYSCORE) dycore diagcore atmcore utilities

core_reg:
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $(PHYSICS) Registry.xml > Registry_processed.xml

core_input_gen:
	if [ ! -e default_inputs ]; then  mkdir default_inputs; fi
	@# Since both the streams.atmosphere file and new stream_list.atmoshpere 
	@# files are regenerated on every make invocation, delete all such existing
	@# files in default_inputs in order ensure that any stale stream_list.atmosphere.
	@# files are removed.  Removal of streams.atmosphere is not necessary since
	@# it will always be overwritten, but we do so to maintain the same behavior
	@# for all stream-related files.
	rm -f default_inputs/streams.atmosphere default_inputs/stream_list.atmosphere.*
	( cd default_inputs; $(NL_GEN) ../Registry_processed.xml namelist.atmosphere in_defaults=true )
	( cd default_inputs; $(ST_GEN) ../Registry_processed.xml streams.atmosphere stream_list.atmosphere. listed in_defaults=true)

gen_includes: core_reg
	(if [ ! -d inc ]; then mkdir -p inc; fi) # To generate *.inc files
	(cd inc; $(REG_PARSE) ../Registry_processed.xml $(CPPFLAGS) )

post_build:
	@#
	@# In the following line, ROOT_DIR is the path to the top-level MPAS-Model
	@# directory.  This if-statement creates a default_inputs directory under
	@# this root directory if it doesn't already exist, e.g.
	@#
	@#   /path/to/my/MPAS-Model/default_inputs
	@#
	@# Call this the "main" default_inputs directory.
	@#
	if [ ! -e $(ROOT_DIR)/default_inputs ]; then mkdir $(ROOT_DIR)/default_inputs; fi
	@# 
	@# Remove streams.atmosphere and any stream_list.atmosphere. files from
	@# ROOT_DIR/default_inputs before copying into this directory the newly-
	@# generated files from the current core's default_inputs directory.  This
	@# ensures that stale streams.atmosphere files from previous make invocations
	@# are removed from ROOT_DIR/default_inputs.
	@#
	rm -f $(ROOT_DIR)/default_inputs/streams.$(CORE) $(ROOT_DIR)/default_inputs/stream_list.$(CORE).*
	@#
	@# The following line copies everything from the default_inputs subidirectory
	@# under the current core's source directory into the main default_inputs
	@# directory.  Since there is a Makefile like this under each core, we call
	@# the former the (current) core's default_inputs directory.  For example, if
	@# this Makefile is for the atmosphere core, it will be located under
	@#
	@#   /path/to/my/MPAS-Model/src/core_atmosphere
	@#
	@# Then the following command will copy everything from
	@#
	@#   /path/to/my/MPAS-Model/src/core_atmosphere/default_inputs
	@#
	@# into
	@#
	@#   /path/to/my/MPAS-Model/default_inputs
	@#
	@# Note that the core's default_inputs directory is guaranteed to not
	@# contain any stale stream_list.atmosphere. files from previous make
	@# invocations (see the rm command in the core_input_gen target).  
	@#
	cp default_inputs/* $(ROOT_DIR)/default_inputs/.
	@#
	@################################################################
	@# The following line will copy each file in the main default_inputs
	@# directory one level up into the top-level MPAS-Model directory, i.e.
	@# ROOT_DIR, but ONLY IF that file doesn't already exist in ROOT_DIR.
	@#
	@# CAUTION:
	@# If files from an old build are already present in ROOT_DIR and a
	@# new build causes them to be updated in ROOT_DIR/default_inputs,
	@# the following line WILL NOT update those files in ROOT_DIR.
	@#
	@# This seems to be intentional behavior, but it's not clear why.
	@#
	@#( cd $(ROOT_DIR)/default_inputs; for FILE in `ls -1`; do if [ ! -e ../$$FILE ]; then cp $$FILE ../.; fi; done )
	@#
	@################################################################
	@#
	@# Update line so it ALWAYS copies over all files from ROOT_DIR/default_inputs
	@# into ROOT_DIR.
	@#
	@#( cd $(ROOT_DIR)/default_inputs; for FILE in `ls -1`; do cp $$FILE ../.; done )
	@#
	@################################################################
	@#
	@# Even better is to make backups of any existing files in ROOT_DIR
	@# (by renaming them as .old001, .old002, etc) and then copying the
	@# file in ROOT_DIR/default_inputs over to ROOT_DIR.
	@#
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

dycore: mpas_atm_dimensions.o $(PHYSCORE)
	( cd dynamics; $(MAKE) all PHYSICS="$(PHYSICS)" )

diagcore: $(PHYSCORE) dycore
	( cd diagnostics; $(MAKE) all PHYSICS="$(PHYSICS)" )

utilities: $(PHYSCORE)
	( cd utils; $(MAKE) all PHYSICS="$(PHYSICS)" )

atmcore: $(PHYSCORE) dycore diagcore $(OBJS)
	ar -ru libdycore.a $(OBJS) dynamics/*.o $(PHYS_OBJS) diagnostics/*.o

mpas_atm_core_interface.o: mpas_atm_core.o

mpas_atm_core.o: dycore diagcore mpas_atm_threading.o mpas_atm_halos.o

mpas_atm_dimensions.o:

clean:
	( cd physics; $(MAKE) clean )
	( cd dynamics; $(MAKE) clean )
	( cd stochastic_physics; $(MAKE) clean )
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
	$(CPP) $(CPPFLAGS) $(PHYSICS) $(CPPINCLUDES) -I./inc $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework -I../operators -I./physics -I./dynamics -I./diagnostics -I./physics/physics_wrf -I./physics/physics_mmm -I./physics/physics_noaa/UGWP -I./stochastic_physics -I../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(PHYSICS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I./inc -I../framework -I../operators -I./physics -I./dynamics -I./diagnostics -I./physics/physics_wrf -I./physics/physics_mmm -I./physics/physics_noaa/UGWP -I./stochastic_physics -I../external/esmf_time_f90
endif
