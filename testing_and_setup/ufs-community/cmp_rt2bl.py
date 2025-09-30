#!/usr/bin/env python

##############################################################################
# Dependencies
##############################################################################
import os
import sys
from os.path import exists
import argparse

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
##############################################################################
def compare_files(dir_bl, dir_rt, filein):
    file_bl = dir_bl+'/'+filein
    file_rt = dir_rt+'/'+filein
    error_count   = 0
    message = 'Comparing ' + file_bl + ' to ' + file_rt
    if (os.path.isfile(file_rt)):
        com = 'nccmp -d ' + file_bl + ' ' + file_rt + ' > logfile.txt'
        result = os.system(com)
        if (result != 0):
            message = message + '  NOT IDENTICAL'
            error_count = error_count + 1
        else:
            message = message + '  PASS'
        # End if
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

    return error_count
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

    # MPAS baseline files
    file_bl_hist = get_files(dir_bl,'history.','.nc')
    file_bl_diag = get_files(dir_bl,'diag.','.nc')

    # Compare MPAS baselines to feature branch
    print('-'*50)
    error_count = 0
    for file_hist in file_bl_hist:
        error_count = error_count + compare_files(dir_bl, dir_rt, file_hist)
    # end for
    for file_diag in file_bl_diag:
        error_count = error_count + compare_files(dir_bl, dir_rt, file_diag)
    # end for
        
    if error_count == 0:
        print("ALL TESTS PASSED, OUTPUT IS IDENTICAL.")
    else:
        print("ALL TESTS PASSED, BUT OUTPUT DIFFERS FROM BASELINE.")
    # end if
    return error_count
# end def

##############################################################################
#
##############################################################################
if __name__ == '__main__':
    error_count = main()
    print(error_count)
# end main
