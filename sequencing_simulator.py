#!/usr/bin/env python3
# Kayla Schimke, Christopher Vollmers

import sys
import os
import time
import glob
import numpy as np
import argparse

def commandline():
    parser = argparse.ArgumentParser(description='Simulates an ONT sequencing run based on the metadata from input files',
                                     add_help=True,prefix_chars='-',
                                     usage='python3 sequencing_simulator.py -i [fast5s] -o [new directory]')
    parser.add_argument('--input','-i', type=str, required=True, help='Directory containing fast5 files.')
    parser.add_argument('--output','-o', type=str, default=os.getcwd(), help='Path to directory to copy fast5s. Defaults to current directory.')

    args = parser.parse_args()
    return args

def main():
    args = commandline()

    if not os.path.exists(args.output):
        os.mkdir(args.output)

    file_list = glob.glob(args.input + '/*.fast5')
    file_list.sort(key=lambda x:os.path.getmtime(x))

    previous = os.path.getmtime(file_list[0])
    avg_time = []
    for f5_path in file_list:
        current = os.path.getmtime(f5_path)
        wait_time = current-previous
        fast5_name = f5_path.split('/')[-1]

        print('Waiting %f seconds to move file %s' %(wait_time,fast5_name))
        time.sleep(wait_time)
        os.system('scp %s %s' %(f5_path,args.output+fast5_name))

        if wait_time > 0:
            avg_time.append(wait_time)

        previous = current

    print('Average time between files: %f seconds' %(np.median(avg_time)))

if __name__ == '__main__':
    main()
