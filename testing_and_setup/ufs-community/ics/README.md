# make_ics: init_atmosphere tests

`.github/workflows/make_ics.yml` builds `init_atmosphere_model` from the branch under test and creates the static, initial-condition, lateral-boundary (and, for GFS, surface-update) files for three 120 km CONUS cases:

| case | source | date | static land use |
|---|---|---|---|
| `conus_gfs` | GFS | 2023-03-10 15Z | 20-class MODIS |
| `conus_raps` | RAP | 2024-08-15 18Z | 20-class MODIS |
| `conus_raps_lake` | RAP | 2024-08-15 18Z | 21-class MODIS with lakes (`config_landuse_data = 'MODIFIED_IGBP_MODIS_NOAH_LAKES'`) |

`conus_gfs` and `conus_raps` share the 20-class static file. The summer RAP date is used because its inputs are small (1.3 GB of ungrib files, against 16 GB for winter RAP). `cases/` also defines `conus_rapw` and `conus_rapw_lake`, which can be run by hand or added to the workflow matrix.

Each case is checked with `check_ics.py`:
- the expected files exist;
- `nlcat` is 21 for `conus_raps_lake` and 20 for the others;
- `islake_lu` is 21 for `conus_raps_lake` and -1 for the others;
- in `conus_raps_lake` the lake points (category 21) exist and all have `landmask` 0; the others have none;
- the main init fields contain no NaN or Inf;
- the LBC file has the expected number of times (17 for GFS over 48 h, 2 for RAP), and GFS has an `sfc_update` file.

The files are uploaded as `ics-<case>` artifacts.

The existing `run_mpas*` workflows still use the staged files from `data/get_data.sh` and are unchanged.

## Layout

```
cases/<case>.sh              case settings (sourced by create_case.sh)
case_files/                  init_atmosphere namelists and streams per step
  step1_static/              20-class static file
  step1_static_lakes/        21-class static file with lakes
  step2_init/, step3_lbc/, step4_sst/   per source and date
scripts/create_case.sh       runs the steps for one case
scripts/fetch_inputs.sh      downloads and unpacks input bundles from a release
check_ics.py                 sanity checks on one case's files
copy_inputs_release.sh       copies the input release from mpas-ic-factory
```

`scripts/` and `case_files/` come from [dustinswales/mpas-ic-factory](https://github.com/dustinswales/mpas-ic-factory), with two changes. `create_case.sh` reads `cases/<case>.sh` instead of `cases.yaml`, so it needs no Python. `fetch_inputs.sh` streams bundle parts into `tar` one at a time, which keeps disk use down.

## Inputs

The workflow downloads its inputs from the release named by `INPUTS_TAG` in `make_ics.yml` (`v8.4.2-1.19`) of the repository it runs in. Copy them once from the factory release, then do the same for any other repository that runs the workflow, such as ufs-community/MPAS-Model. By default only the bundles the workflow uses are copied (grids, aerosols, geographic data including the lakes dataset, summer RAP, GFS and its SST: about 14 GB). Bundles already in the release are skipped, so re-running adds only what is missing. Add `-a` for every bundle (about 30 GB), or `-b "..."` to choose:

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
