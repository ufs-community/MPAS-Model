#!/usr/bin/env bash
#
# Create MPAS static / init / lbc / sfc_update files for one case in cases/.
# Copied from dustinswales/mpas-ic-factory (scripts/create_case.sh); case settings
# come from cases/<CASE>.sh instead of cases.yaml.
# Portable version of barlage/mpas_testcase create_case/create_mpas_testcase_caseNN.sh
# (GSL code and namelists):
# no Slurm, no site paths. Runs on a laptop, a GitHub runner or an HPC node.
#
# Usage:
#   scripts/create_case.sh -c CASE -e INIT_ATMOSPHERE_EXE -i INPUT_DIR -o OUTPUT_DIR [options]
#
# Options:
#   -c CASE        case id: a file cases/CASE.sh (e.g. conus_gfs)
#   -e EXE         init_atmosphere_model executable
#   -i INPUT_DIR   input tree (see scripts/fetch_inputs.sh):
#                    grids/<mesh>.grid.nc, grids/<mesh>.graph.info
#                    geog/mpas_static/...                      (static step only)
#                    ungrib/<source>.<yyyymmddhh>/<PREFIX>:<date>
#                    ungrib/sst.<yyyymmddhh>/SST:<date>         (gfs only)
#                    aerosols/QNWFA_QNIFA_SIGMA_MONTHLY.dat     (use_climo_aerosols)
#   -o OUTPUT_DIR  final files go to OUTPUT_DIR/<case_dir>/case_files
#   -n NTASKS      MPI tasks (default: nproc)
#   -s STEPS       comma list from static,init,lbc,sst (default: all that apply)
#   -S STATIC_DIR  directory holding an existing mpas.*.static.nc (and ugwp_oro_data)
#                  for this case; used instead of running the static step
#   -w WORK_DIR    run directories (default: OUTPUT_DIR/<case_dir>)
#   -k             keep the step run directories (default: remove after success)
#
# Environment:
#   MPIRUN   launcher prefix, task count appended (default "mpirun -n"; e.g. "srun -n")
#   GPMETIS  gpmetis executable (default "gpmetis")
#   DESCRIBE_EXTRA  file whose contents are appended to case_files/description
#
set -euo pipefail

script_home="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

die() { echo "ERROR: $*" >&2; exit 1; }
banner() {
  echo
  echo "################################################################"
  echo "# $*"
  echo "################################################################"
}

case_id="" exe="" input_dir="" output_dir="" ntasks="" steps="" static_dir="" work_root="" keep="false"
while getopts "c:e:i:o:n:s:S:w:kh" opt; do
  case $opt in
    c) case_id=$OPTARG ;;
    e) exe=$OPTARG ;;
    i) input_dir=$OPTARG ;;
    o) output_dir=$OPTARG ;;
    n) ntasks=$OPTARG ;;
    s) steps=$OPTARG ;;
    S) static_dir=$OPTARG ;;
    w) work_root=$OPTARG ;;
    k) keep="true" ;;
    h) sed -n '2,32p' "$0"; exit 0 ;;
    *) die "bad option (see -h)" ;;
  esac
done
[ -n "$case_id" ]    || die "-c CASE is required"
[ -n "$exe" ]        || die "-e EXE is required"
[ -n "$input_dir" ]  || die "-i INPUT_DIR is required"
[ -n "$output_dir" ] || die "-o OUTPUT_DIR is required"
[ -x "$exe" ]        || die "executable not found: $exe"
exe="$(cd "$(dirname "$exe")" && pwd)/$(basename "$exe")"
input_dir="$(cd "$input_dir" && pwd)"
mkdir -p "$output_dir"; output_dir="$(cd "$output_dir" && pwd)"
ntasks=${ntasks:-$(nproc)}
MPIRUN=${MPIRUN:-"mpirun -n"}
GPMETIS=${GPMETIS:-gpmetis}

# Case settings: case_id domain resolution source season date mesh case_dir
#                static_key use_climo_aerosols
domain="" resolution="" source="" date="" mesh="" case_dir="" use_climo_aerosols=""
static_variant="" static_nml_dir="step1_static"
case_file="$script_home/cases/$case_id.sh"
[ -f "$case_file" ] || die "unknown case '$case_id' (no $case_file)"
# shellcheck source=/dev/null
. "$case_file"

yyyy=${date:0:4} mm=${date:4:2} dd=${date:6:2} hh=${date:8:2}
datestring="${yyyy}-${mm}-${dd}_${hh}.00.00"
prefix="mpas.$domain.$resolution"
static_file="$prefix.static.nc"
ugwp_file="$prefix.ugwp_oro_data.nc"
init_file="$prefix.$source.init.$datestring.nc"
lbc_file="$prefix.$source.lbc.$datestring.nc"
sst_file="$prefix.$source.sfc_update.$datestring.nc"

has_lbc=false; [ "$domain" = "conus" ] && has_lbc=true
has_sst=false; [ "$source" = "gfs" ] && has_sst=true

if [ -z "$steps" ]; then
  steps="init"
  [ -z "$static_dir" ] && steps="static,$steps"
  $has_lbc && steps="$steps,lbc"
  $has_sst && steps="$steps,sst"
fi
want() { [[ ",$steps," == *",$1,"* ]]; }
want lbc && ! $has_lbc && die "lbc step only applies to conus (case $case_id is $domain)"
want sst && ! $has_sst && die "sst step only applies to gfs source (case $case_id is $source)"
if want init && ! want static && [ -z "$static_dir" ]; then
  die "init without static: pass -S STATIC_DIR or include the static step"
fi

case_root="$output_dir/$case_dir"
work=${work_root:-$case_root}
final="$case_root/case_files"
mkdir -p "$work" "$final"
work="$(cd "$work" && pwd)"

grid_file="$input_dir/grids/$mesh.grid.nc"
graph_file="$input_dir/grids/$mesh.graph.info"
ungrib_met="$input_dir/ungrib/$source.$date"
ungrib_sst="$input_dir/ungrib/sst.$date"
geog_dir="$input_dir/geog/mpas_static"
aerosol_file="$input_dir/aerosols/QNWFA_QNIFA_SIGMA_MONTHLY.dat"

# Check every input up front so a missing bundle fails before any long step runs.
if want static; then
  [ -s "$grid_file" ] || die "missing $grid_file"
  [ -d "$geog_dir" ]  || die "missing $geog_dir"
fi
if want init || want lbc; then
  [ -d "$ungrib_met" ] || die "missing $ungrib_met"
  if [ "$use_climo_aerosols" = "true" ]; then
    [ -s "$aerosol_file" ] || die "missing $aerosol_file (use_climo_aerosols=true)"
  fi
fi
if want sst; then [ -d "$ungrib_sst" ] || die "missing $ungrib_sst"; fi

# MPAS needs a large stack (the MPAS user guide says "ulimit -s unlimited");
# the 8 MB default segfaults init_atmosphere on some compilers/inputs.
ulimit -s unlimited 2>/dev/null || ulimit -s "$(ulimit -H -s)" 2>/dev/null || true
export OMP_STACKSIZE=${OMP_STACKSIZE:-512M}

banner "Configuration settings"
cat <<EOF
case:               $case_id ($case_dir)
steps:              $steps
executable:         $exe
domain/resolution:  $domain/$resolution (mesh $mesh)
static namelist:    case_files/$static_nml_dir${static_variant:+  (variant $static_variant)}
input source:       $source  date $datestring
input dir:          $input_dir
work dir:           $work
output dir:         $final
ntasks:             $ntasks  (launcher: $MPIRUN)
stack limit:        $(ulimit -s)
use_climo_aerosols: $use_climo_aerosols
EOF

# ------------------------------------------------------------------
# helpers
# ------------------------------------------------------------------
partition() {  # link graph.info.part.N into the current directory
  [ "$ntasks" -gt 1 ] || return 0
  local part="$work/partitions/$mesh.graph.info.part.$ntasks"
  if [ ! -s "$part" ]; then
    if [ -s "$input_dir/grids/$mesh.graph.info.part.$ntasks" ]; then
      part="$input_dir/grids/$mesh.graph.info.part.$ntasks"
    else
      [ -s "$graph_file" ] || die "missing $graph_file (needed to partition for $ntasks tasks)"
      mkdir -p "$work/partitions"
      cp "$graph_file" "$work/partitions/$mesh.graph.info"
      (cd "$work/partitions" && $GPMETIS -minconn -contig -niter=200 "$mesh.graph.info" "$ntasks" >/dev/null) \
        || die "gpmetis failed for $mesh with $ntasks parts"
    fi
  fi
  ln -sf "$part" "graph.info.part.$ntasks"
}

stage_namelist() {  # $1 = step directory name, $2 = optional dated subdirectory
  local src="$script_home/case_files/$1${2:+/$2}"
  [ -d "$src" ] || die "no namelist directory $src"
  cp "$src"/namelist.init_atmosphere "$src"/streams.init_atmosphere .
}

run_init_atmosphere() {  # $1 = expected output file
  ln -sf "$exe" init_atmosphere_model
  rm -f log.init_atmosphere.*
  local t0=$SECONDS rc=0
  $MPIRUN "$ntasks" ./init_atmosphere_model || rc=$?
  echo "init_atmosphere_model finished in $((SECONDS - t0)) s (rc=$rc)"
  local errs=(log.init_atmosphere.*.err)
  [ -e "${errs[0]}" ] || errs=()
  if [ $rc -ne 0 ] || [ ! -s "$1" ]; then
    # Full .out log first, then the warnings/errors, so the cause ends up at the
    # bottom of the CI log. ("SEAICE_FRACTIONAL/OMLD file not found" lines are
    # informational: those are optional inputs.)
    if [ -e log.init_atmosphere.0000.out ]; then
      echo "---- tail -80 log.init_atmosphere.0000.out"; tail -80 log.init_atmosphere.0000.out
    fi
    for f in "${errs[@]}"; do echo "---- $f"; cat "$f"; done
    echo "---- ERROR/CRITICAL lines"
    grep -hE "ERROR|CRITICAL" log.init_atmosphere.* 2>/dev/null | sort -u || echo "(none in logs)"
    [ -s "$1" ] || die "$1 was not created in $PWD (rc=$rc)"
    die "init_atmosphere_model exited with $rc in $PWD"
  elif [ ${#errs[@]} -gt 0 ]; then
    echo "---- warnings (run succeeded)"; cat "${errs[@]}"
  fi
}

save_log() {  # $1 = step run dir, $2 = name in case_files
  [ -e "$1/log.init_atmosphere.0000.out" ] && cp "$1/log.init_atmosphere.0000.out" "$final/$2"
  return 0
}

link_aerosols() {
  if [ "$use_climo_aerosols" = "true" ]; then
    [ -s "$aerosol_file" ] || die "missing $aerosol_file (use_climo_aerosols=true)"
    ln -sf "$aerosol_file" .
  fi
}

banner "Start processing time: $(date -u)"

# ------------------------------------------------------------------
# step 1: static
# ------------------------------------------------------------------
static_run="$work/step1_static"
if want static; then
  banner "Start creating: $static_file"
  [ -s "$grid_file" ] || die "missing $grid_file"
  [ -d "$geog_dir" ]  || die "missing $geog_dir"
  rm -rf "$static_run"; mkdir -p "$static_run"; cd "$static_run"
  stage_namelist "$static_nml_dir"
  ln -sf "$geog_dir" mpas_static
  ln -sf "$grid_file" grid.nc
  partition
  run_init_atmosphere static.nc
  [ -s ugwp_oro_data.nc ] || die "ugwp_oro_data.nc was not created"
  cp static.nc "$final/$static_file"
  cp ugwp_oro_data.nc "$final/$ugwp_file"
  save_log "$static_run" log.static
fi

# Static file for the later steps.
if want init; then
  if [ -n "$static_dir" ]; then
    static_src="$static_dir/$static_file"
    [ -s "$static_src" ] || die "no $static_file in $static_dir"
    cp "$static_src" "$final/$static_file"
    [ -s "$static_dir/$ugwp_file" ] || die "no $ugwp_file in $static_dir"
    cp "$static_dir/$ugwp_file" "$final/$ugwp_file"
    [ -e "$static_dir/log.static" ] && cp "$static_dir/log.static" "$final/log.static"
  fi
  [ -s "$final/$static_file" ] || die "no static file $final/$static_file"
fi

# ------------------------------------------------------------------
# step 2: init
# ------------------------------------------------------------------
init_run="$work/step2_init"
if want init; then
  banner "Start creating: $init_file"
  [ -d "$ungrib_met" ] || die "missing $ungrib_met"
  rm -rf "$init_run"; mkdir -p "$init_run"; cd "$init_run"
  stage_namelist step2_init "$source.$date"
  partition
  shopt -s nullglob
  met0=("$ungrib_met"/*"${yyyy}-${mm}-${dd}_${hh}")
  shopt -u nullglob
  [ ${#met0[@]} -gt 0 ] || die "no intermediate file for ${yyyy}-${mm}-${dd}_${hh} in $ungrib_met"
  ln -sf "${met0[@]}" .
  ln -sf "$final/$static_file" static.nc
  link_aerosols
  run_init_atmosphere init.nc
  mv init.nc "$final/$init_file"
  save_log "$init_run" log.init
fi

# ------------------------------------------------------------------
# step 3: lbc (conus)
# ------------------------------------------------------------------
lbc_run="$work/step3_lbc"
if want lbc; then
  banner "Start creating: $lbc_file"
  [ -s "$final/$init_file" ] || die "lbc needs $final/$init_file"
  rm -rf "$lbc_run"; mkdir -p "$lbc_run"; cd "$lbc_run"
  stage_namelist step3_lbc "$source.$date"
  partition
  ln -sf "$ungrib_met"/* .
  ln -sf "$final/$init_file" init.nc
  link_aerosols
  run_init_atmosphere lbc.nc
  mv lbc.nc "$final/$lbc_file"
  save_log "$lbc_run" log.lbc
fi

# ------------------------------------------------------------------
# step 4: sfc_update (gfs)
# ------------------------------------------------------------------
sst_run="$work/step4_sst"
if want sst; then
  banner "Start creating: $sst_file"
  [ -s "$final/$init_file" ] || die "sst needs $final/$init_file"
  [ -d "$ungrib_sst" ] || die "missing $ungrib_sst"
  rm -rf "$sst_run"; mkdir -p "$sst_run"; cd "$sst_run"
  stage_namelist step4_sst "$source.$date"
  partition
  ln -sf "$ungrib_sst"/SST* .
  ln -sf "$final/$init_file" init.nc
  run_init_atmosphere sfc_update.nc
  mv sfc_update.nc "$final/$sst_file"
  save_log "$sst_run" log.sst
fi

# ------------------------------------------------------------------
# provenance
# ------------------------------------------------------------------
cd "$final"
{
  echo "These files were created with:"
  echo "case:               $case_id"
  echo "created (UTC):      $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "host:               $(hostname)"
  echo "executable:         $exe"
  echo "executable sha256:  $(sha256sum "$exe" | cut -d' ' -f1)"
  echo "launcher:           $MPIRUN $ntasks"
  echo "use_climo_aerosols: $use_climo_aerosols"
  echo "input dir:          $input_dir"
  [ -s "$input_dir/INPUTS_TAG" ] && echo "inputs tag:         $(cat "$input_dir/INPUTS_TAG")"
  [ -n "${DESCRIBE_EXTRA:-}" ] && [ -s "$DESCRIBE_EXTRA" ] && cat "$DESCRIBE_EXTRA"
  true
} > description
sha256sum mpas.*.nc > SHA256SUMS

if [ "$keep" = "false" ]; then
  rm -rf "$static_run" "$init_run" "$lbc_run" "$sst_run" "$work/partitions"
fi

banner "Successful completion of case: $case_dir ($steps)"
ls -l "$final"
banner "End processing time: $(date -u)"
