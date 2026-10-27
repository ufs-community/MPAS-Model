# make_ics: init_atmosphere tests

`.github/workflows/make_ics.yml` builds `init_atmosphere_model` from the branch under test and creates the static, initial-condition, lateral-boundary and surface-update files for four 120 km CONUS cases:

| case | source | date | static land use |
|---|---|---|---|
| `conus_gfs` | GFS | 2023-03-10 15Z | 20-class MODIS |
| `conus_raps` | RAP | 2024-08-15 18Z | 20-class MODIS |
| `conus_rapw` | RAP | 2024-02-02 18Z | 20-class MODIS |
| `conus_rapw_lake` | RAP | 2024-02-02 18Z | 21-class MODIS with lakes (`config_landuse_data = 'MODIFIED_IGBP_MODIS_NOAH_LAKES'`) |

Each case is checked with `check_ics.py`:
- the expected files exist;
- `nlcat` is 21 for `conus_rapw_lake` and 20 for the others;
- `islake_lu` is 21 for `conus_rapw_lake` and -1 for the others;
- in `conus_rapw_lake` the lake points (category 21) exist and all have `landmask` 0; no other case has any;
- the main init fields contain no NaN or Inf;
- the LBC file has the expected number of times.

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

The workflow downloads its inputs from the `init-inputs-v1.0.0` release of the repository it runs in. Copy them once from the factory release, then do the same for any other repository that runs the workflow, such as ufs-community/MPAS-Model:

```bash
bash testing_and_setup/ufs-community/ics/copy_inputs_release.sh dustinswales/MPAS-Model /path/to/work
```

## Running by hand

```bash
bash scripts/fetch_inputs.sh -r dustinswales/MPAS-Model -t init-inputs-v1.0.0 -d inputs \
    grids geog aerosols ungrib-rap.2024020218
MPIRUN="mpirun -np" bash scripts/create_case.sh -c conus_rapw_lake \
    -e /path/to/init_atmosphere_model -i inputs -o out -n 4
python check_ics.py out/conus.120km.rap.2024020218.lakes/case_files --source rap --lakes --lbc-times 2
```
