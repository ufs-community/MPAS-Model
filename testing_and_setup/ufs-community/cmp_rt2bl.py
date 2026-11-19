#!/usr/bin/env python

##############################################################################
# Dependencies
##############################################################################
import os
import sys
import shutil
import subprocess
from os.path import exists
import argparse

# nccmp exit statuses: 0 = identical, 1 = files differ; anything else (2 = fatal error, e.g. a
# file that cannot be read; 126/127 = nccmp cannot be run) is an error, not a difference.
NCCMP_DIFFER = 1
# Exit status of this script when a comparison could not be made. Differences in the output
# (and missing feature-branch files) are reported but do not change the exit status.
EXIT_ERROR = 2

##############################################################################
# Command line arguments
##############################################################################
parser = argparse.ArgumentParser()
parser.add_argument('-drt', '--dir_rt', help='MPAS run directory for testing',   required=True)
parser.add_argument('-dbl', '--dir_bl', help='MPAS run directory for baselines', required=True)

def parse_args():
    args   = parser.parse_args()
    dir_rt = args.dir_rt
    dir_bl = args.dir_bl
    return (dir_rt, dir_bl)
# end def

##############################################################################
# Procedure to return <file_list> from a given <directory> provided the
# file <prefix> and <suffix>.
##############################################################################
def get_files(directory, prefix, suffix):
    file_list = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(suffix):
                if file.startswith(prefix):
                    file_list.append(file)
                # end if
            # end if
        # end for
    # end for
    return file_list
# end def

##############################################################################
# Procedure to compare <filein> in <dir_bl> to <dir_rt>.
# NOTE.This procedure assumes the baseline <filein> exists.
# Returns (number of files that differ or are missing, number of comparisons that failed).
##############################################################################
def compare_files(dir_bl, dir_rt, filein):
    file_bl = dir_bl+'/'+filein
    file_rt = dir_rt+'/'+filein
    error_count   = 0
    failed_count  = 0
    message = 'Comparing ' + file_bl + ' to ' + file_rt
    if (os.path.isfile(file_rt)):
        try:
            result = subprocess.run(['nccmp', '-d', file_bl, file_rt],
                                    capture_output=True, text=True)
            status = result.returncode
            output = (result.stdout + result.stderr).strip()
        except OSError as e:
            status = None
            output = str(e)
        # end try
        if (status == 0):
            message = message + '  PASS'
        elif (status == NCCMP_DIFFER):
            message = message + '  NOT IDENTICAL'
            error_count = error_count + 1
            # Show the first differences nccmp found
            lines = output.splitlines()
            message = message + ''.join('\n    ' + l for l in lines[:10])
            if len(lines) > 10:
                message = message + '\n    ... (' + str(len(lines) - 10) + ' more)'
            # end if
        else:
            message = message + '  ERROR: nccmp ' + \
                      ('could not be run' if status is None else 'exited with status ' + str(status))
            message = message + ''.join('\n    ' + l for l in output.splitlines()[:20])
            failed_count = failed_count + 1
        # end if
    else:
        if not exists(file_rt):
            message = message + '  MISSING testing file:  ' + file_rt
        # end if
        if not exists(file_bl):
            message = message + '  MISSING baseline file: ' + file_bl
        # end if
        error_count = error_count + 1
    # end if
    print(message)

    return error_count, failed_count
# end def

##############################################################################
# Main program
# Given two MPAS run directories, <dir_bl> and <dir_rt>, compare all of the
# history and diagnostic files within
# <dir_bl> contains MPAS output from the unmodified authoratative codebase,
# the "baseline".
# <dir_rt> contains MPAS output from the "feature" branch being propsed to the
# authoratative gsl/develop for inclusion.
#
##############################################################################
def main():
    # Get command line arguments.
    (dir_rt, dir_bl) = parse_args()

    if shutil.which('nccmp') is None:
        print('ERROR: nccmp not found on PATH; cannot compare output.')
        return EXIT_ERROR
    # end if

    # MPAS baseline files
    file_bl_hist = get_files(dir_bl,'history.','.nc')
    file_bl_diag = get_files(dir_bl,'diag.','.nc')

    # Compare MPAS baselines to feature branch
    print('-'*50)
    error_count  = 0
    failed_count = 0
    for filein in file_bl_hist + file_bl_diag:
        (differ, failed) = compare_files(dir_bl, dir_rt, filein)
        error_count  = error_count  + differ
        failed_count = failed_count + failed
    # end for

    if failed_count > 0:
        print("COMPARISON FAILED: " + str(failed_count) + " file(s) could not be compared.")
        return EXIT_ERROR
    elif error_count == 0:
        print("ALL TESTS PASSED, OUTPUT IS IDENTICAL.")
    else:
        print("ALL TESTS PASSED, BUT OUTPUT DIFFERS FROM BASELINE.")
    # end if
    return 0
# end def

##############################################################################
#
##############################################################################
if __name__ == '__main__':
    sys.exit(main())
# end main
