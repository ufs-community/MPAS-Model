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
    no_plots   = args.no_plots
    return (dir_rt, dir_bl, no_plots)

#
def main():
    #
    (dir_rt, dir_bl, no_plots) = parse_args()

    #
    error_count = 0
    for run in run_list:
        print(run_list)
        for root, dirs, files in os.walk(dir_rt):
            for file in files:
                if file.endswith('.nc'):
                    print(file)
                # end if
            # end for
        # end for
    # end for

#
if __name__ == '__main__':
    main()
