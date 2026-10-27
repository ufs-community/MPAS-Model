#!/usr/bin/env bash
#
# Copy the make_ics input bundles from the mpas-ic-factory release into a release
# of this repository (the workflow reads release INPUTS_TAG of its own repository).
#
# Usage:
#   bash copy_inputs_release.sh [-s SRC_REPO] [-t SRC_TAG] [-T DEST_TAG] [-b "BUNDLES" | -a] DEST_REPO WORK_DIR
#     -s   source repository   (default dustinswales/mpas-ic-factory)
#     -t   source release tag  (default v1.0.0)
#     -T   destination tag     (default v8.4.2-1.19, as INPUTS_TAG in make_ics.yml)
#     -b   bundles to copy (default: what make_ics.yml uses, about 14 GB:
#          "grids aerosols geog ungrib-rap.2024081518 ungrib-gfs.2023031015
#          ungrib-sst.2023031015"; geog includes geog-* bundles). Bundles already in
#          the destination release are skipped, so -b can add to it.

#     -a   copy every bundle in the source release (about 30 GB)
#   e.g.
#     bash copy_inputs_release.sh dustinswales/MPAS-Model /scratch/$USER/copy
#
# Needs gh with write access to DEST_REPO (gh auth login, or GH_TOKEN). Assets are
# copied one at a time (download, upload, delete), so WORK_DIR needs room for the
# largest single asset (< 2 GB) only. Already-copied assets with the same checksum are
# skipped, so the script can be re-run after an interruption.
[ -n "${BASH_VERSION:-}" ] || { echo "ERROR: run with bash: bash $0 ..." >&2; exit 1; }
set -euo pipefail

die() { echo "ERROR: $*" >&2; exit 1; }
usage="usage: $0 [-s SRC_REPO] [-t SRC_TAG] [-T DEST_TAG] [-b \"BUNDLES\" | -a] DEST_REPO WORK_DIR"

src_repo="dustinswales/mpas-ic-factory" src_tag="v1.0.0" dest_tag="v8.4.2-1.19"
bundles="grids aerosols geog ungrib-rap.2024081518 ungrib-gfs.2023031015 ungrib-sst.2023031015"
all=false
while getopts "s:t:T:b:a" opt; do
  case $opt in
    s) src_repo=$OPTARG ;;
    t) src_tag=$OPTARG ;;
    T) dest_tag=$OPTARG ;;
    b) bundles=$OPTARG ;;
    a) all=true ;;
    *) die "$usage" ;;
  esac
done
shift $((OPTIND - 1))
[ $# -eq 2 ] || die "$usage"
dest_repo=$1 work=$2

command -v gh >/dev/null || die "gh CLI not found"
gh auth status >/dev/null 2>&1 || die "gh is not logged in (set GH_TOKEN or run gh auth login)"
mkdir -p "$work"; work="$(cd "$work" && pwd)"
cd "$work"

gh release download "$src_tag" -R "$src_repo" -p MANIFEST.sha256 -O src_all.sha256 --clobber \
  || die "cannot read MANIFEST.sha256 from $src_repo@$src_tag"

# Select the assets of the requested bundles (BUNDLE.tar.gz or BUNDLE.tar.gz.part-xx;
# "geog" also selects every geog-<dataset> bundle).
if $all; then
  cp src_all.sha256 src.sha256
else
  : > src.sha256
  for b in $bundles; do
    pat="^${b//./\\.}\\.tar\\.gz(\\.part-[a-z]+)?\$"
    [ "$b" = geog ] && pat="^geog(-[^ ]+)?\\.tar\\.gz(\\.part-[a-z]+)?\$"
    n0=$(wc -l < src.sha256)
    PAT="$pat" awk '$2 ~ ENVIRON["PAT"]' src_all.sha256 >> src.sha256
    [ "$(wc -l < src.sha256)" -gt "$n0" ] || die "bundle '$b' not in $src_repo@$src_tag"
  done
fi
echo "Copying $(wc -l < src.sha256) assets from $src_repo@$src_tag to $dest_repo@$dest_tag"

if gh release view "$dest_tag" -R "$dest_repo" >/dev/null 2>&1; then
  gh release download "$dest_tag" -R "$dest_repo" -p MANIFEST.sha256 -O dest.sha256 --clobber 2>/dev/null \
    || : > dest.sha256
else
  gh release create "$dest_tag" -R "$dest_repo" --title "init_atmosphere test inputs $dest_tag" \
    --notes "Inputs for .github/workflows/make_ics.yml, copied from $src_repo@$src_tag."
  : > dest.sha256
fi

n=$(wc -l < src.sha256); i=0
while read -r sum asset; do
  i=$((i + 1))
  if grep -qxF "$sum  $asset" dest.sha256; then
    echo "[$i/$n] $asset: already in $dest_repo@$dest_tag"
    continue
  fi
  echo "[$i/$n] $asset"
  gh release download "$src_tag" -R "$src_repo" -p "$asset" --clobber
  echo "$sum  $asset" | sha256sum -c --quiet - || die "$asset: checksum mismatch after download"
  gh release upload "$dest_tag" -R "$dest_repo" --clobber "$asset"
  rm -f "$asset"
  # record progress so a re-run skips it
  awk -v f="$asset" '$2 != f' dest.sha256 > dest.tmp; echo "$sum  $asset" >> dest.tmp; mv dest.tmp dest.sha256
done < src.sha256

# Merged manifest: every source entry, plus any destination entries not in the source.
{ awk 'NR == FNR { s[$2] = 1; next } !($2 in s)' src.sha256 dest.sha256; cat src.sha256; } > MANIFEST.sha256
gh release upload "$dest_tag" -R "$dest_repo" --clobber MANIFEST.sha256
echo "Done: $dest_repo@$dest_tag has $(wc -l < MANIFEST.sha256) assets listed."
