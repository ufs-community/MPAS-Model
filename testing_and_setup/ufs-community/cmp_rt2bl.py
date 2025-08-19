#!/usr/bin/env python

##############################################################################
#
# This script compares MPAS RT output to baselines.
#
##############################################################################
import os
import sys
from rt_test_cases import run_list
from os.path import exists
import argparse
#from plot_scm_out import plot_results

#
parser = argparse.ArgumentParser()
parser.add_argument('-drt',  '--dir_rt',   help='Directory containing SCM RT output',              required=True)
parser.add_argument('-dbl',  '--dir_bl',   help='Directory containing SCM RT baselines',           required=True)
parser.add_argument('-np',   '--no_plots', help='flag to turn off generation of difference plots', required=False, action='store_true')

#
def parse_args():
    args     = parser.parse_args()
    dir_rt   = args.dir_rt
    dir_bl   = args.dir_bl
    no_plots = args.no_plots
    return (dir_rt, dir_bl, no_plots)

##############################################################################
# Procedure to return <file_list> in a provided <directory> for the
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
    error_message = ''
    if (os.path.isfile(file_rt)):
        com = 'nccmp -d ' + file_bl + ' ' + file_rt + ' > logfile.txt'
        print('Comparing ',file_bl,' to ',file_rt)
        result = os.system(com)
        if (result != 0):
            message = "  NOT IDENTICAL"
            error_count = error_count + 1
        else:
            message = "  PASS"
        # End if
        print(message)
        # end if
        ierr = 0
    else:
        if not exists(file_rt):
            message = "  MISSING testing file: " + file_rt
        # end if
        if not exists(file_bl):
            message = "  MISSING baseline file: " + file_bl
        # end if
        error_count = error_count + 1
    # end if

    return error_count
# end def

def main():
    #
    (dir_rt, dir_bl, no_plots) = parse_args()

    #
    for run in run_list:
        print(run_list)

        # MPAS baseline files
        file_bl_hist = get_files(dir_bl,'history.','.nc')
        file_bl_diag = get_files(dir_bl,'diag.','.nc')

        # Compare baselines to regression_test.
        print('-'*50)
        for file_hist in file_bl_hist:
            ierr = compare_files(dir_bl, dir_rt, file_list)
        # end for
        for file_diag in file_bl_diag:
            ierr = ierr + compare_files(dir_bl, dir_rt, file_diag)
        # end for
        
        if ierr == 0:
            print("ALL TESTS PASSED, OUTPUT IS IDENTICAL.")
        else:
            print("ALL TESTS PASSED, BUT OUTPUT DIFFERS FROM BASELINE.")
        # end if
#
if __name__ == '__main__':
    main()
