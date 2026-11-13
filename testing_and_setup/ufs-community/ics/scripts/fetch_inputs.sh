#!/usr/bin/env bash
#
# Download input bundles from a GitHub release and unpack them into INPUT_DIR.
# Copied from dustinswales/mpas-ic-factory (scripts/fetch_inputs.sh), with parts
# streamed into tar one at a time to keep disk use down on CI runners.
#
# Usage:
#   scripts/fetch_inputs.sh -r OWNER/REPO -t TAG -d INPUT_DIR BUNDLE [BUNDLE ...]
#
# BUNDLE names (assets are BUNDLE.tar.gz, or BUNDLE.tar.gz.part-aa, -ab, ... when split):
#   grids                 -> INPUT_DIR/grids/
#   geog                  -> INPUT_DIR/geog/mpas_static/
#                            (also fetches every geog-<dataset> bundle in the release,
#                             see tools/add_geog_dataset.sh)
#   aerosols              -> INPUT_DIR/aerosols/
#   ungrib-<src>.<date>   -> INPUT_DIR/ungrib/<src>.<date>/   (e.g. ungrib-gfs.2023031015)
#
# Every asset is checked against MANIFEST.sha256 from the same release.
# Uses the gh CLI when available (GH_TOKEN), otherwise plain https downloads
# (public repositories only).
set -euo pipefail

die() { echo "ERROR: $*" >&2; exit 1; }

repo="" tag="" input_dir=""
while getopts "r:t:d:" opt; do
  case $opt in
    r) repo=$OPTARG ;;
    t) tag=$OPTARG ;;
    d) input_dir=$OPTARG ;;
    *) die "usage: $0 -r OWNER/REPO -t TAG -d INPUT_DIR BUNDLE..." ;;
  esac
done
shift $((OPTIND - 1))
[ -n "$repo" ] && [ -n "$tag" ] && [ -n "$input_dir" ] || die "-r, -t and -d are required"
[ $# -gt 0 ] || die "no bundles requested"

mkdir -p "$input_dir"
input_dir="$(cd "$input_dir" && pwd)"
dl="$input_dir/.download"
mkdir -p "$dl"

use_gh=false
if command -v gh >/dev/null 2>&1 && [ -n "${GH_TOKEN:-${GITHUB_TOKEN:-}}" ]; then use_gh=true; fi

fetch() {  # fetch asset names into $dl
  if $use_gh; then
    local args=()
    for a in "$@"; do args+=(-p "$a"); done
    gh release download "$tag" -R "$repo" -D "$dl" --clobber "${args[@]}"
  else
    for a in "$@"; do
      curl -fL --retry 5 --retry-delay 10 -sS -o "$dl/$a" \
        "https://github.com/$repo/releases/download/$tag/$a"
    done
  fi
}

if [ -e "$input_dir/INPUTS_TAG" ] && [ "$(cat "$input_dir/INPUTS_TAG")" != "$repo@$tag" ]; then
    die "$input_dir contains inputs from another release; use a fresh input directory"
fi
echo "$repo@$tag" > "$input_dir/INPUTS_TAG"
 
[ -s "$dl/MANIFEST.sha256" ] || fetch MANIFEST.sha256

# "geog" also brings in extra datasets added later as geog-<dataset> bundles.
bundles=("$@")
for b in "$@"; do
  if [ "$b" = geog ]; then
    mapfile -t extra < <(awk '{print $2}' "$dl/MANIFEST.sha256" \
                         | sed -nE 's/^(geog-[^ ]+)\.tar\.gz(\.part-[a-z]+)?$/\1/p' | sort -u)
    bundles+=("${extra[@]}")
  fi
done

for bundle in "${bundles[@]}"; do
  marker="$input_dir/.done.$bundle"
  if [ -e "$marker" ]; then echo "$bundle: already unpacked"; continue; fi

  mapfile -t assets < <(awk '{print $2}' "$dl/MANIFEST.sha256" \
                         | grep -E "^${bundle//./\\.}\.tar\.gz(\.part-[a-z]+)?$" | sort)
  [ ${#assets[@]} -gt 0 ] || die "bundle '$bundle' not in release $repo@$tag"

  # One part at a time: download, verify, feed to tar, delete. Peak disk use is the
  # unpacked data plus one part, instead of every part plus the unpacked data.
  echo "$bundle: downloading and unpacking ${#assets[@]} asset(s)"
  for a in "${assets[@]}"; do
    fetch "$a" >&2
    (cd "$dl" && awk -v f="$a" '$2 == f' MANIFEST.sha256 | sha256sum -c --quiet -) >&2 \
      || { echo "ERROR: $a: checksum mismatch" >&2; exit 1; }
    cat "$dl/$a"
    rm -f "$dl/$a"
  done | tar -xzf - -C "$input_dir" || die "$bundle: download or unpack failed"
  touch "$marker"
done

echo "$repo@$tag" > "$input_dir/INPUTS_TAG"
du -sh "$input_dir"/* 2>/dev/null || true
