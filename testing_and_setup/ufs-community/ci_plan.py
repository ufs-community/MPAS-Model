#!/usr/bin/env python3
"""Make the build and run matrices for .github/workflows/mpas_ci.yml from ci_tests.json.

Writes `builds=<json>` and `runs=<json>` to $GITHUB_OUTPUT (or prints them), each a matrix
of the form {"include": [...]}. Only the builds some test needs are listed.

If DOCKER_TREE is set (the git tree hash of docker/, first 12 characters), each toolchain
uses the image <image>-<DOCKER_TREE> that bld_mpas_images.yaml builds from those Dockerfiles,
and it is an error if that image does not exist in the registry. Otherwise the images in
ci_tests.json are used as they are.

Writes image_<toolchain>=<image> for each toolchain to $GITHUB_OUTPUT as well.

Usage: ci_plan.py [ci_tests.json]
"""
import json
import os
import subprocess
import sys

here = os.path.dirname(os.path.abspath(__file__))
config_file = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "ci_tests.json")
with open(config_file) as f:
    config = json.load(f)
toolchains, builds, tests = config["toolchains"], config["builds"], config["tests"]


def image_exists(image):
    for _ in range(3):  # retry, so a registry hiccup is not reported as a missing image
        try:
            if subprocess.run(["docker", "manifest", "inspect", image], capture_output=True,
                              timeout=60).returncode == 0:
                return True
        except (OSError, subprocess.TimeoutExpired):
            pass
    return False


tree = os.environ.get("DOCKER_TREE", "")
if tree:
    missing = []
    for tc, cfg in toolchains.items():
        cfg["image"] = f"{cfg['image']}-{tree}"
        print(f"{tc}: {cfg['image']}")
        if not image_exists(cfg["image"]):
            missing.append(cfg["image"])
    if missing:
        msg = (f"No CI image built from this docker/ (tree {tree}): {', '.join(missing)}. "
               "Run the 'MPAS-A Base Images Build' workflow (bld_mpas_images.yaml) on a "
               "branch with these Dockerfiles, then re-run this workflow.")
        print(f"::error::{msg}")
        sys.exit(1)

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
lines += [f"image_{tc}={cfg['image']}" for tc, cfg in toolchains.items()]
if os.environ.get("GITHUB_OUTPUT"):
    with open(os.environ["GITHUB_OUTPUT"], "a") as f:
        f.write("\n".join(lines) + "\n")
else:
    print("\n".join(lines))
print(f"{len(build_list)} builds, {len(runs)} runs")
for b in build_list:
    print(f"  build  {b['toolchain']:7s} {b['variant']}")
for r in runs:
    print(f"  run    {r['toolchain']:7s} {r['test']}")
