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

def get_files(dirIN,prefix):
    file_list = []
    for root, dirs, files in os.walk(dirIN):
        for file in files:
            if file.endswith('.nc'):
                if file.startswith(prefix):
                    file_list.append(file)
                # end if
            # end if
        # end for
    # end for
    return file_list
# end def

def compare_files(dir_bl, dir_rt, files):
    file_bl = dir_bl+'/'+files
    file_rt = dir_rt+'/'+files
    if (os.path.isfile(file_rt)):
        com = 'nccmp -d ' + file_bl + ' ' + file_rt + ' > logfile.txt'
        print('Comparing ',file_bl,' to ',file_rt)
        result = os.system(com)
        if (result != 0):
            message = "  Output DIFFERS from baseline."
        else:
            message = "  Output IDENTICAL to baseline."
        # End if
        print(message)
        # end if
        ierr = 0
    else:
        print("ERROR: Cannot find file for comparision, ",file_rt)
        ierr = 1
    # end if
    return ierr
# end def

def main():
    #
    (dir_rt, dir_bl, no_plots) = parse_args()

    #
    for run in run_list:
        print(run_list)

        # MPAS history files (baselines)
        file_bl_hist = get_files(dir_bl,'history.')

        # MAPS diagnsotic files (baselines)
        file_bl_diag = get_files(dir_bl,'diag.')

        # Compare baselines to regression_test.
        for files in file_bl_hist:
            ierr = compare_files(dir_bl, dir_rt, files)
        # end for

        for files in file_bl_diag:
            ierr = compare_files(dir_bl, dir_rt, files)
        # end for
#
if __name__ == '__main__':
    main()
