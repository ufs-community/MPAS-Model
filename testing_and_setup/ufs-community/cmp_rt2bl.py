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
    args    = parser.parse_args()
    dir_rt  = args.dir_rt 
    dir_bl  = args.dir_bl
    no_plots   = args.no_plots
    return (dir_rt, dir_bl, no_plots)

#
def main():
    #
    (dir_rt, dir_bl, no_plots) = parse_args()

    #
    error_count = 0
    for run in run_list:
        file_rt = dir_rt + "/" + run["case"]+"_"+run["suite"]+"/output.nc"
        file_bl = dir_bl + "/" + run["case"]+"_"+run["suite"]+"/output.nc"
        if exists(file_rt) and exists(file_bl):
            com = "cmp "+file_rt+" "+file_bl+" > logfile.txt"
            result = os.system(com)
            if (result != 0):
                message = "Output for "+run["case"]+"_"+run["suite"]+ " DIFFERS from baseline."
                if (not no_plots):
                    message += " Difference plots will be created."
                print(message)
                error_count = error_count + 1
            else:
                print("Output for "+run["case"]+"_"+run["suite"]+ " is IDENTICAL to baseline")
            # end if
        else:
            if not exists(file_rt):
                print("Output for "+run["case"]+"_"+run["suite"]+ " is MISSING from output")
            # end if
            if not exists(file_bl):
                print("Output for "+run["case"]+"_"+run["suite"]+ " is MISSING from baseline")
            # end if
            error_count = error_count + 1
        # end if
    # end for
