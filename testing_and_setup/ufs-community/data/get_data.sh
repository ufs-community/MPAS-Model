#!/bin/bash

# Function to display help message
print_help() {
    echo "get_data.sh: testing_and_setup/ufs-community/data/get_data.sh [-v,--verbose]"
    echo "    Script for downloading/extracting the Physics lookup tables and MPAS ICs."
    echo ""
    echo "Options:"
    echo "    -v, --verbose    Turn on wget verbose output."
    echo "    --help           Show this help message and exit."
}

verbose="-q"
# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --help)
            print_help
            exit 0
            ;;
        -v|--verbose)
            verbose="-v"
            ;;
        *)
            echo "Unknown option: $1"
            print_help
            exit 1
            ;;
    esac
    shift
done

set -ex

if [[ $(uname -s) == Darwin ]]; then
  if [[ $(sw_vers -productVersion) < 12.3 ]]; then
    MYDIR=$(cd "$(dirname "$(greadlink -f -n "${BASH_SOURCE[0]}" )" )" && pwd -P)
  else
    MYDIR=$(cd "$(dirname "$(readlink -f -n "${BASH_SOURCE[0]}" )" )" && pwd -P)
  fi
else
  MYDIR=$(cd "$(dirname "$(readlink -f -n "${BASH_SOURCE[0]}" )" )" && pwd -P)
fi
BASEDIR=$MYDIR/../../..


# Change to directory containing theinput data, download and extract archive
cd $BASEDIR/testing_and_setup/ufs-community/data/

# Get TEMPO data
wget ${verbose} https://github.com/ufs-community/MPAS-Model/releases/download/MPAS-v8.3.1-2.14/tempo_data.tar
mkdir -p tables/tempo/
mv tempo_data.tar tables/tempo/
cd tables/tempo
tar -xvf tempo_data.tar
rm tempo_data.tar
cd ../../

# Get Thompson data
wget ${verbose} https://github.com/ufs-community/MPAS-Model/releases/download/MPAS-v8.3.1-2.13/thompson_data.tar
mkdir -p tables/thompson/
mv thompson_data.tar tables/thompson/
cd tables/thompson
tar -xvf thompson_data.tar
rm thompson_data.tar
cd ../../

# Get UGW data
wget ${verbose} https://github.com/ufs-community/MPAS-Model/releases/download/MPAS-v8.3.1-2.13/ugw_data.tar
mkdir -p tables/ugw/
mv ugw_data.tar tables/ugw/
cd tables/ugw
tar -xvf ugw_data.tar
rm ugw_data.tar
cd ../../

# Get MPAS case data
wget ${verbose} https://github.com/ufs-community/MPAS-Model/releases/download/MPAS-v8.3.1-2.13/mpas_data.tar
mkdir -p ics
mv mpas_data.tar ics/
cd ics
tar -xvf mpas_data.tar
rm mpas_data.tar
cd ../
