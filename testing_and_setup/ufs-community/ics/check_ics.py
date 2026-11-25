#!/usr/bin/env python3
"""Sanity checks on the files made by scripts/create_case.sh for one make_ics case.

  check_ics.py CASE_FILES_DIR --source {gfs,rap} [--lakes] [--noahmp] [--soil-levels N]
               [--lbc-times N]

Checks:
  - static, ugwp_oro_data, init and lbc files exist (and sfc_update for gfs)
  - nlcat is 21 with --lakes, 20 otherwise
  - islake_lu is 21 with --lakes, -1 (no lake category) otherwise
  - with --lakes: lake points (ivgtyp == 21) exist and are all water (landmask 0);
    without: no point has category 21
  - the init file has --soil-levels soil levels (default 4 with --noahmp, 9 otherwise)
  - with --noahmp: the Noah-MP static fields (soilcomp, soilcl1-4) are in the static
    and init files and hold valid values over land, and every soil level of tslb
    and smois holds a plausible value over land (150-350 K, 0-1)
  - no NaN/Inf in the main init fields
  - the lbc file has --lbc-times time records (when given)
Exit status 1 if any check fails.
"""
import argparse
import glob
import os
import sys

import numpy as np
from netCDF4 import Dataset

NOAHMP_FIELDS = ["soilcomp", "soilcl1", "soilcl2", "soilcl3", "soilcl4"]
INIT_FIELDS = ["theta", "rho", "u", "w", "qv", "surface_pressure", "skintemp", "tslb", "smois"]


def one(pattern, what, failures):
    files = sorted(glob.glob(pattern))
    if len(files) != 1:
        failures.append(f"expected one {what} file ({os.path.basename(pattern)}), found {len(files)}")
        return None
    return files[0]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("case_files")
    p.add_argument("--source", choices=["gfs", "rap"], required=True)
    p.add_argument("--lakes", action="store_true")
    p.add_argument("--noahmp", action="store_true")
    p.add_argument("--soil-levels", type=int)
    p.add_argument("--lbc-times", type=int)
    a = p.parse_args()

    d = a.case_files
    failures = []
    ok = []

    static = one(f"{d}/mpas.*.static.nc", "static", failures)
    one(f"{d}/mpas.*.ugwp_oro_data.nc", "ugwp_oro_data", failures)
    init = one(f"{d}/mpas.*.init.*.nc", "init", failures)
    lbc = one(f"{d}/mpas.*.lbc.*.nc", "lbc", failures)
    if a.source == "gfs":
        one(f"{d}/mpas.*.sfc_update.*.nc", "sfc_update", failures)

    want_nlcat = 21 if a.lakes else 20
    want_islake = 21 if a.lakes else -1

    for path in [x for x in (static, init) if x]:
        name = os.path.basename(path)
        with Dataset(path) as nc:
            nlcat = len(nc.dimensions["nlcat"]) if "nlcat" in nc.dimensions else None
            if nlcat != want_nlcat:
                failures.append(f"{name}: nlcat = {nlcat}, expected {want_nlcat}")
            else:
                ok.append(f"{name}: nlcat = {nlcat}")
            if "islake_lu" in nc.variables:
                islake = int(np.asarray(nc.variables["islake_lu"][...]).ravel()[0])
                if islake != want_islake:
                    failures.append(f"{name}: islake_lu = {islake}, expected {want_islake}")
                else:
                    ok.append(f"{name}: islake_lu = {islake}")
            else:
                failures.append(f"{name}: no islake_lu variable")

            lu = np.asarray(nc.variables["ivgtyp"][:]).ravel()
            lm = np.asarray(nc.variables["landmask"][:]).ravel()
            nlake = int((lu == 21).sum())
            if a.lakes:
                if nlake == 0:
                    failures.append(f"{name}: no lake points (ivgtyp == 21)")
                elif np.any(lm[lu == 21] != 0):
                    failures.append(f"{name}: {int((lm[lu == 21] != 0).sum())} of {nlake} lake points have landmask != 0")
                else:
                    ok.append(f"{name}: {nlake} lake points, all landmask = 0")
            elif nlake:
                failures.append(f"{name}: {nlake} points with ivgtyp == 21 without the lakes option")

            if a.noahmp:
                land = lm == 1
                for v in NOAHMP_FIELDS:
                    if v not in nc.variables:
                        failures.append(f"{name}: no {v} variable (Noah-MP static field)")
                        continue
                    var = nc.variables[v]
                    x = np.asarray(var[:], dtype=np.float64)
                    # land points only: move the nCells axis first, then select
                    x = np.moveaxis(x, var.dimensions.index("nCells"), 0)
                    xl = x[land]
                    if not np.all(np.isfinite(xl)):
                        failures.append(f"{name}: {v} has NaN/Inf over land")
                    elif v == "soilcomp" and not np.any((xl > 0) & (xl <= 100)):
                        failures.append(f"{name}: soilcomp has no values in (0, 100] over land")
                    elif v != "soilcomp" and not np.any(xl >= 1):
                        failures.append(f"{name}: {v} has no soil texture class >= 1 over land")
                    else:
                        ok.append(f"{name}: {v} valid over land")

    if init:
        with Dataset(init) as nc:
            want_soil = a.soil_levels or (4 if a.noahmp else 9)
            nsoil = len(nc.dimensions["nSoilLevels"]) if "nSoilLevels" in nc.dimensions else None
            if nsoil != want_soil:
                failures.append(f"{os.path.basename(init)}: nSoilLevels = {nsoil}, expected {want_soil}")
            else:
                ok.append(f"{os.path.basename(init)}: nSoilLevels = {nsoil}")
            if a.noahmp:
                land = np.asarray(nc.variables["landmask"][:]).ravel() == 1
                for v, lo, hi in (("tslb", 150.0, 350.0), ("smois", 0.0, 1.0)):
                    if v not in nc.variables:
                        continue
                    var = nc.variables[v]
                    x = np.moveaxis(np.asarray(var[:], dtype=np.float64), var.dimensions.index("nCells"), 0)[land]
                    x = x.reshape(x.shape[0], -1)
                    bad = [k + 1 for k in range(x.shape[1]) if not np.all((x[:, k] >= lo) & (x[:, k] <= hi))]
                    if bad:
                        failures.append(f"{os.path.basename(init)}: {v} outside [{lo:g}, {hi:g}] over land at soil level(s) {bad}")
                    else:
                        ok.append(f"{os.path.basename(init)}: {v} in [{lo:g}, {hi:g}] over land at all {x.shape[1]} soil levels")
            for v in INIT_FIELDS:
                if v not in nc.variables:
                    failures.append(f"{os.path.basename(init)}: missing required init field {v}")
                    continue
                x = np.asarray(nc.variables[v][:], dtype=np.float64)
                bad = int((~np.isfinite(x)).sum())
                if bad:
                    failures.append(f"{os.path.basename(init)}: {v} has {bad} NaN/Inf values")
            ok.append(f"{os.path.basename(init)}: no NaN/Inf in {', '.join(v for v in INIT_FIELDS if v in nc.variables)}")

    if lbc and a.lbc_times is not None:
        with Dataset(lbc) as nc:
            nt = len(nc.dimensions["Time"])
        if nt != a.lbc_times:
            failures.append(f"{os.path.basename(lbc)}: {nt} time records, expected {a.lbc_times}")
        else:
            ok.append(f"{os.path.basename(lbc)}: {nt} time records")

    for line in ok:
        print(f"  ok    {line}")
    for line in failures:
        print(f"  FAIL  {line}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
