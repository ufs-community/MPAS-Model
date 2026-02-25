#!/usr/bin/env bash
# shellcheck disable=SC1090,1091
HOMEmpas="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
source "${HOMEmpas}/ush/detect_machine.sh"

COMPILER="${COMPILER:-intel}"
# ==============================================================================
usage() {
  set +x
  echo
  echo "Usage: $0 -j <num> -h"
  echo
  echo "  -j  number of build jobs               DEFAULT: 8"
  echo "  -c  additional CMAKE options"
  echo "  -f  force a clean build"
  echo "  -h  display this message and quit"
  echo
  exit 1
}

while getopts "c:j:hf" opt; do
  case ${opt} in
    c)
      CMAKE_OPTS=${OPTARG}
      ;;
    j)
      BUILD_JOBS=${OPTARG}
      ;;
    f)
      CLEAN_BUILD=YES
      ;;
    h|\?|:)
      usage
      ;;
  esac
done

cd "${HOMEmpas}" || exit 1

module purge                      
module use "${HOMEmpas}/modulefiles"
module load "mpas/${MACHINE}.${COMPILER}"
module list

BUILD_DIR=${BUILD_DIR:-${HOMEmpas}/build}
if [[ ${CLEAN_BUILD} == 'YES' ]]; then
  [[ -d ${BUILD_DIR} ]] && rm -rf ${BUILD_DIR}
elif [[ -d ${BUILD_DIR} ]]; then
  printf "Build directory (${BUILD_DIR}) already exists\n"
  printf "Please choose what to do:\n\n"
  printf "[r]emove the existing directory\n"
  printf "[c]ontinue building in the existing directory\n"
  printf "[q]uit this build script\n"
  read -p "Choose an option (r/c/q):" choice
  case ${choice} in
    [Rr]* ) rm -rf ${BUILD_DIR} ;;
    [Cc]* ) ;;
        * ) exit ;;
  esac
fi
mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR} || exit 1

BUILD_JOBS=${BUILD_JOBS:-8}
CMAKE_OPTS+=' -DMPAS_DOUBLE_PRECISION=OFF'
cmake ${CMAKE_OPTS} -DMPAS_CORES="init_atmosphere;atmosphere" ..
make -j $BUILD_JOBS
