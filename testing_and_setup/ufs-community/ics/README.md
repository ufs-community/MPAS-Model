# Initial conditions for MPAS CI (init_atmosphere tests)

Step 1 of `.github/workflows/mpas_ci.yml` builds `init_atmosphere_model` from the branch under test and creates the static, initial-condition, lateral-boundary (and, for GFS, surface-update) files for six 120 km CONUS cases:

| case | source | date | static land use |
|---|---|---|---|
| `conus_gfs` | GFS | 2023-03-10 15Z | 20-class MODIS |
| `conus_gfs_lake` | GFS | 2023-03-10 15Z | 21-class MODIS with lakes |
| `conus_gfs_noahmp` | GFS | 2023-03-10 15Z | 20-class MODIS; ICs for Noah-MP (see below) |
| `conus_gfs_noahmp10` | GFS | 2023-03-10 15Z | 20-class MODIS; ICs for Noah-MP with 10 soil layers |
| `conus_raps` | RAP | 2024-08-15 18Z | 20-class MODIS |
| `conus_raps_lake` | RAP | 2024-08-15 18Z | 21-class MODIS with lakes (`config_landuse_data = 'MODIFIED_IGBP_MODIS_NOAH_LAKES'`) |

`conus_gfs`, `conus_gfs_noahmp`, `conus_gfs_noahmp10` and `conus_raps` share the 20-class static file, and `conus_gfs_lake` and `conus_raps_lake` share the 21-class one, so there are two static jobs. The summer RAP date is used because its inputs are small (1.3 GB of ungrib files, against 16 GB for winter RAP). `cases/` also defines `conus_rapw` and `conus_rapw_lake`, which can be run by hand or added to the workflow matrix.

`conus_gfs_noahmp` makes ICs for Noah-MP, which needs 4 soil levels (the GSL namelists use 9, for RUC). Its init, LBC and sfc_update namelists are the NCAR ones from barlage/mpas_testcase (55 levels, 30 km model top, `case_files/step*/gfs.2023031015.noahmp`), as used for the `mesoscale_reference_noahmp` run tests. The Noah-MP static fields (`soilcomp`, `soilcl1`-`soilcl4`) are already in every static file, because `config_noahmp_static` defaults to true. A case picks its own init/LBC/sfc_update namelists with `nml_variant` in `cases/<case>.sh`.

`conus_gfs_noahmp10` is the same with `config_nsoillevels = 10`. `init_atmosphere` sets the 10 layer thicknesses (`dzs` = 0.05, 0.10, 0.10, 0.10, 0.30, 0.35, 0.50, 1.00, 1.00, 1.50 m, 5 m in all) in `src/core_atmosphere/physics/mpas_atmphys_initialize_real.F` and interpolates soil temperature and moisture to them as for 4 layers. GFS soil data and the deep soil temperature reach only 3 m, so the deepest layer (centered at 4.25 m) takes the deepest input values.

Each case is checked with `check_ics.py`:
- the expected files exist;
- `nlcat` is 21 for the `_lake` cases and 20 for the others;
- `islake_lu` is 21 for the `_lake` cases and -1 for the others;
- in the `_lake` cases the lake points (category 21) exist and all have `landmask` 0; the others have none;
- the init file has 4 soil levels for `conus_gfs_noahmp`, 10 for `conus_gfs_noahmp10` and 9 for the others;
- for the Noah-MP cases, the static and init files hold valid `soilcomp` and `soilcl1`-`soilcl4` values over land, and every soil level of `tslb` and `smois` is plausible over land;
- the main init fields contain no NaN or Inf;
- the LBC file has the expected number of times (17 for GFS over 48 h, 2 for RAP), and GFS has an `sfc_update` file.

The files are uploaded as `ics-<case>` artifacts.

Step 2 of the workflow (the `run` job) runs `atmosphere_model` from the branch under test and from a baseline on these files: the `mpasdev.*` run tests use `conus_gfs_noahmp` (NCAR levels), `ufscommunity.hrrrv5.gfs.winter` and the stochastic tests use `conus_gfs`, and `ufscommunity.hrrrv5.rap.summer` uses `conus_raps`. Only the physics tables still come from `data/get_data.sh` (`--tables-only`).

## Layout

```
cases/<case>.sh              case settings (sourced by create_case.sh)
case_files/                  init_atmosphere namelists and streams per step
  step1_static/              20-class static file
  step1_static_lakes/        21-class static file with lakes
  step2_init/, step3_lbc/, step4_sst/   per source and date (<source>.<date>[.<nml_variant>])
scripts/create_case.sh       runs the steps for one case
scripts/fetch_inputs.sh      downloads and unpacks input bundles from a release
check_ics.py                 sanity checks on one case's files
copy_inputs_release.sh       copies the input release from mpas-ic-factory
```

`scripts/` and `case_files/` come from [dustinswales/mpas-ic-factory](https://github.com/dustinswales/mpas-ic-factory), with two changes. `create_case.sh` reads `cases/<case>.sh` instead of `cases.yaml`, so it needs no Python. `fetch_inputs.sh` streams bundle parts into `tar` one at a time, which keeps disk use down.

## Inputs

The workflow downloads its inputs from the release named by `INPUTS_TAG` in `mpas_ci.yml` of the repository it runs in. Copy them once from the factory release, then do the same for any other repository that runs the workflow, such as ufs-community/MPAS-Model. By default only the bundles the workflow uses are copied (grids, aerosols, geographic data including the lakes dataset, summer RAP, GFS and its SST: about 14 GB). Bundles already in the release are skipped, so re-running adds only what is missing. Add `-a` for every bundle (about 30 GB), or `-b "..."` to choose:

```bash
bash testing_and_setup/ufs-community/ics/copy_inputs_release.sh dustinswales/MPAS-Model /path/to/work
```

## Running by hand

```bash
bash scripts/fetch_inputs.sh -r dustinswales/MPAS-Model -t v8.4.2-1.19 -d inputs \
    grids geog aerosols ungrib-rap.2024081518
MPIRUN="mpirun -np" bash scripts/create_case.sh -c conus_raps_lake \
    -e /path/to/init_atmosphere_model -i inputs -o out -n 4
python check_ics.py out/conus.120km.rap.2024081518.lakes/case_files --source rap --lakes --lbc-times 2
```
