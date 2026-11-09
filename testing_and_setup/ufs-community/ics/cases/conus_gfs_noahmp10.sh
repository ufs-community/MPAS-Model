# make_ics case conus_gfs_noahmp10 (sourced by scripts/create_case.sh)
# As conus_gfs_noahmp, with 10 soil layers for Noah-MP instead of 4:
#   dzs = 0.05 0.10 0.10 0.10 0.30 0.35 0.50 1.00 1.00 1.50 m (5 m deep), set in
#   init_atmosphere (src/core_atmosphere/physics/mpas_atmphys_initialize_real.F)
# when config_nsoillevels = 10. Layers centered below 3 m (the deepest GFS input
# level, where the deep soil temperature is used) take the deepest input values.
domain=conus
resolution=120km
mesh=conus.120km
source=gfs
date=2023031015
case_dir=conus.120km.gfs.2023031015.noahmp10
static_key=conus-120km
static_variant=
static_nml_dir=step1_static
nml_variant=noahmp10
use_climo_aerosols=true
