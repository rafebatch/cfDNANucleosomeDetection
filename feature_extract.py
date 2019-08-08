import sys
import os
import time
import numpy as np
from optparse import OptionParser, OptionGroup
from itertools import groupby
from operator import itemgetter
from bisect import bisect_left

start_time = time.time()
"""
Program to get number of chromosomes per bin.
Used for analysis/differentiation between healthy/cancerous samples.
Probably best to use a different approach, but the code might be useful.

Sorry for lack of documentation, basically just gets number of nucleosomes
per user defined window along the reference genome and outputs to a text file.
"""

"""
Option handling
"""

parser = OptionParser("%prog [options] bed_file1 bed_file2 bed_file3 ...")
parser.add_option("-w", "--window", dest="window_size", help="Option to set the size of the windows.", default = 1000000, type="int")
parser.add_option("-s", "--status", dest="status", help="Specify whether the sample is healthy or unhealthy", default = None)
parser.add_option("-c", "--cancer", dest="cancer", help="Specify the type of cancer a sample has or if healthy, specify none.", default=None)

(options, args) = parser.parse_args()
if options.status == None:
    sys.stderr.write("Please specify whether the input file is a healthy sample\n")
    sys.exit()
if options.cancer == None:
    sys.stderr.write("Please specify the cancer type.\n")
    sys.exit()

"""
File handling
"""

files = []
name_list = []
fnames = []

file_list = args
for each in file_list:
    files.append(each)
try:
    for each in files:
        name_list.append(each.split('.')[1])
        fnames.append(each.split('.')[0])
except:
    sys.stderr.write("Please specify a file in the format: filename.bed\n")
    sys.exit()
for each in name_list:
    if each == "bed":
        continue
    else:
        sys.stderr.write("Please include files in the correct format (bed).\n")
        sys.exit()

def normalize(tot_list):
    norm_vals = (tot_list - np.mean(tot_list)) / np.std(tot_list)
    norm_list = norm_vals.tolist()
    return norm_list


def find_closest(sublist, win_list, window_size, end_pos):
    window_pos = window_size

    while window_pos < end_pos:
        
        pos = bisect_left(sublist, window_pos)

        win_list.append(sublist[:pos])
        sublist = sublist[pos:]
    
        window_pos += window_size        

"""
General functionality
"""

window_size = options.window_size

tot_list = []
window_list = []
file_features = []
min_max_list = []
for each in files:
    with open(each, "r") as f:
        sublist = []
        count = 0
        for line in f:
            vals = line.split() 
            sublist.append(int(vals[2]))
            count += 1
        
        min_max_list.append(int(vals[2]))
        tot_list.append(sublist)

min_max = max(min_max_list)


for each in tot_list:
    window_list = []
    feature_list = []
    find_closest(each, window_list, window_size, min_max)
    for group in window_list:
        feature_list.append(len(group))
    file_features.append(feature_list)

print("num_files=%d\twindow_size=%d\tstatus=%s\tcancer=%s" % (len(args), window_size, options.status, options.cancer))
for sub in file_features:
     for each in sub:
          print(str(each))
