#!/usr/bin/env python3
"""Make the build and run matrices for .github/workflows/mpas_ci.yml from ci_tests.json.

Writes `builds=<json>` and `runs=<json>` to $GITHUB_OUTPUT (or prints them), each a matrix
of the form {"include": [...]}. Only the builds some test needs are listed.

Usage: ci_plan.py [ci_tests.json]
"""
import json
import os
import sys

here = os.path.dirname(os.path.abspath(__file__))
config_file = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "ci_tests.json")
with open(config_file) as f:
    config = json.load(f)
toolchains, builds, tests = config["toolchains"], config["builds"], config["tests"]

errors = []
runs, needed = [], {}
for t in tests:
    for key in ("bl", "ft"):
        if t[key] not in builds:
            errors.append(f"{t['test']}: unknown build {t[key]!r}")
    if not os.path.isdir(os.path.join(here, "cases", t["case"])):
        errors.append(f"{t['test']}: no run config cases/{t['case']}")
    for tc in t["toolchains"]:
        if tc not in toolchains:
            errors.append(f"{t['test']}: unknown toolchain {tc!r}")
            continue
        runs.append({"toolchain": tc, "image": toolchains[tc]["image"],
                     "mpirun": toolchains[tc]["mpirun"],
                     **{k: v for k, v in t.items() if k != "toolchains"}})
        for b in (t["bl"], t["ft"]):
            needed[(tc, b)] = True
if errors:
    sys.exit("ci_tests.json: " + "; ".join(errors))

build_list = [{"toolchain": tc, "variant": b, "image": toolchains[tc]["image"],
               "make_target": toolchains[tc]["make_target"],
               "cpu_flag": toolchains[tc].get("cpu_flag", ""), **builds[b]}
              for (tc, b) in needed if b in builds]

out = {"builds": {"include": build_list}, "runs": {"include": runs}}
lines = [f"{k}={json.dumps(v, separators=(',', ':'))}" for k, v in out.items()]
if os.environ.get("GITHUB_OUTPUT"):
    with open(os.environ["GITHUB_OUTPUT"], "a") as f:
        f.write("\n".join(lines) + "\n")
print(f"{len(build_list)} builds, {len(runs)} runs")
for b in build_list:
    print(f"  build  {b['toolchain']:7s} {b['variant']}")
for r in runs:
    print(f"  run    {r['toolchain']:7s} {r['test']}")
