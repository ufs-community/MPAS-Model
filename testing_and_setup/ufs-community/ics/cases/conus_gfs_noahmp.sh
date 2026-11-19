# make_ics case conus_gfs_noahmp (sourced by scripts/create_case.sh)
# GFS ICs for Noah-MP: 4 soil levels (Noah-MP needs 4; the GSL namelists use 9 for RUC).
# The init/lbc/sfc_update namelists are the NCAR ones from barlage/mpas_testcase
# (55 levels, 30 km model top), as used for mesoscale_reference_noahmp. The static
# file is the 20-class one shared with conus_gfs; it already holds the Noah-MP fields
# (soilcomp, soilcl1-4), since config_noahmp_static defaults to true.
domain=conus
resolution=120km
mesh=conus.120km
source=gfs
date=2023031015
case_dir=conus.120km.gfs.2023031015.noahmp
static_key=conus-120km
static_variant=
static_nml_dir=step1_static
nml_variant=noahmp
use_climo_aerosols=true
