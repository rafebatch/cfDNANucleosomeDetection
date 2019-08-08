import sys
import os
import statistics as st
import scipy.signal
import numpy as np
from bx.intervals.intersection import Intersecter, Interval


"""
sorry if you had to install some stuff

IMPORTANT FOR RUNNING: to save output results, specify output file in command line with  > results.bed
e.g. $ python peak_call.py example.wig > example.bed
"""

"""
File handling
"""

## Ensure user specifies a wig file to be read
if len(sys.argv) < 2:
    sys.stderr.write("Please specificy the wig file to be read from.\n")
    sys.exit()
else:
    infile_name = sys.argv[-1]
    ## Determine if file ends in ".wig"
    try:
        check = infile_name.split('.')[1]
    except:
        sys.stderr.write("Please specify a file in the format: filename.wig\n")
        sys.exit()
    ## if wig file specified, open file
    if check == "wig":
        infile = open(sys.argv[-1], "r")
    else:
        sys.stderr.write("Please include a file in the correct format (wig).\n")
        sys.exit()
"""
Methods
"""

def window_maker(total_list, window_size):
    """
    Yield successive window-sized chunks from the total_list.
    """
    for i in range(0, len(total_list), window_size):
        yield total_list[i:i + window_size]

def set_median(total_list):
    """
    Normalizes WPS score in each sublist created by window_maker.
    The sublist is normalized to a median of choice -- if median is set to median,
    The sublist has a median of 0. By default it is set to be normalized by median * .90
    to add some confidence to WPS scores (i.e. higher scores present, middle range scores less so)
    """
    i = 0
    for sublist in total_list:
        median = st.median(sublist)
        ## set median * whatever you want to add/subtract from amount added to normalize.
        median = median * .90
        if median >= 0:
            sublist[:] = [x - median for x in sublist]
        else:
            sublist[:] = [x + abs(median) for x in sublist]

        total_list[i] = sublist
        i += 1

def find_contig(interval, orig_list, index_start, pot_list):
    """
    Method to find the largest contiguous sum within each nucleosome interval.
    This determines the area within the interval of highest scores (i.e. where the most protection occurs)
    The initial position of this region, the middle of this region, the end of this region,
    the max, and the max position are recorded for each potential nucleosome in that order:
    [start, middle, end, max, max position]

    This results in a 2D array of nucleosome info (pot_list)
    """
    if (interval.start - interval.end) % 2 != 0:
        sublist = orig_list[interval.start - index_start:interval.end - index_start]
    else:
        sublist = orig_list[interval.start - index_start:interval.end - index_start+1]

    # index of median in original list
    initial_pos = interval.start

    """ Finds msub-array with maximum sum and returns the maximum sum and sub-array -- i.e max contig subarray in each region"""
    max_start, max_stop = 0, 1  # Start, Stop of maximum sub-array
    curr = 0                    # pointer to current array
    max_sum = orig_list[0]      # Sum of maximum array
    current_sum = 0             # sum of current array

    for i, elem in enumerate(sublist):
        current_sum += elem

        if max_sum < current_sum:
            max_sum = current_sum 
            max_start = curr
            max_stop = i  + 1
        if current_sum < 0:
            current_sum = 0
            curr = i + 1

    ## get max val, find where it is in sublist
    max_val = max(sublist)
    max_pos, = np.where(sublist == max_val)

    ## this list is taken as a parameter and updated for each potential nucleosome.
    pot_list.append([max_start + initial_pos, initial_pos + (max_stop + max_start) // 2, initial_pos + max_stop, max_val, initial_pos + max_pos[0]])

def find_min_region(contig_list, orig_list, start_pos, chrom_num):
    """
    This method actually goes unused -- couldn't get great results,
    but could lead to more confident results if tweaks made.

    --The last stage of the process, this method introduces some metrics to
    add to the confidence of the nucleosome position data.

    This method finds the minimum score between each nucleosome. So for each nucleosome, a minimum is recorded for the min
    to the left of the nucleosome, between the current nucleosome and the previous nucleosome, and to the right of the nucleosome,
    between the current and next nucleosome.

    The length of this region is used to determine the probability that a nucleosome can exist in this space.
    If the length of this region is large enough (i.e. greater than span_size, which is 200 bp by default), and if
    the max of the potential nucleosome interval is 45 bp away from both minimums, than accept the interval as
    a nucleosome -- the start/end points are printed to the output bed file.
    """
    
    orig_list = list(orig_list)
    subl = orig_list[:contig_list[0][4] - start_pos]
    min_1 = subl.index(min(subl))
    
    ## get the position of the minimum minimum of the first sublist
    min_1_pos = orig_list.index(min(orig_list[:contig_list[0][4] - start_pos]))
    
    for i, read in enumerate(contig_list):

        index = read[4] - start_pos

        ## if there is another nucleosome, get the minimum between current nucl/next nucl
        if i < len(contig_list) - 1:
            next_index = contig_list[i+1][4] - start_pos
        
        ## if not, get last minimum available; if more nucleosomes present
        ## move to next nucleosome after recording minimums
        if i == len(contig_list) - 1:
            min_2_pos = orig_list.index(min(orig_list[index:]))
        else:
            subl = orig_list[index:next_index]
            min_2_pos = subl.index(min(subl)) + index

        ## use mins surrounding interval to get total span of region
        len_win = abs(min_2_pos - min_1_pos)
       
        ## not actually sure why I did this lol
        vals = [min_1_pos, min_2_pos]
        
        ## parameter to change for varying confidence results
        span_size = 200

        ## minimum position between each read should span AT LEAST 200 bp if a nucleosome is within range
        if index + 45 < vals[1] and index - 45 > vals[0] and len_win > span_size:
            sys.stdout.write("chr" + str(chrom_num) + "\t" + str(read[0]) + "\t" + str(read[2]) + "\n")

        min_1_pos = min_2_pos

def calc_nucl_distance(contig_list, smooth_list, start_pos, chrom_num):
    """
    The REAL last stage in the process. This method performs an alternative way of boosting confidence.
    The position of the maximum protection score for each potential nucleosome region is retrieved (max_index). 
    Two indices (90 base pairs upstream of the max value and 90 base pairs downstream of the max value) are
    retrieved. If the normalized score is below 0 in both of these positions (90 base pairs behind and ahead),
    then we can say that this max falls within a window 180 bp or less.

    The confidence can probably be boosted using the middle position of the nucleosome interval and setting a 
    minimum paramter.
    """
    smooth_list = list(smooth_list)
    for read in contig_list:
        max_index = read[4] - start_pos
        if smooth_list[max_index + 90] < 0 and smooth_list[max_index - 90] < 0:
            sys.stdout.write("chr" + str(chrom_num) + "\t" + str(read[0]) + "\t" + str(read[2]) + "\n")


"""
General Functionality
"""
## add functionality for multiple chromosomes -- not done

## initial position of read
pos = 0

## initialize start_base -- do I even need to do this in python
start_base = 0

## initialize end_base as well for some reason
end_base = 0

## initialize this to false, think thats default but oh well
started = False

## create the Intersecter -- list that holds intervals (start, end)
nucl = Intersecter()

## more initialization
count = 0
thresh = 0

## get all WPS scores from the wig file in order
values = infile.readlines()


## set the window size for median normalization -- this can be played with for varying results
window_size = 1000

## get rid of the first line of information and process that info
info_str = values.pop(0)

## first line grants us which chromosome is being read, as well as the starting pos of the reads
chrom_num = int(info_str.split()[1].split('=')[1].split('r')[1])
startpos = int(info_str.split()[2].split('=')[1])

## convert all read scores to ints
values = list(map(int, values))

## create the subwindows
windows = list(window_maker(values, window_size))

## normalize the list
set_median(windows)

## recreate a complete, undivided list full of normalized values per 1000 bp window
scaled_list = [j for sub in windows for j in sub]

## implement scipy's savitsky golay filter
SGF_window_size = 21
SGF_order = 2

## use a savgol filter to smooth the resulting list -- suggested by OG paper
smooth_list = scipy.signal.savgol_filter(scaled_list, SGF_window_size, SGF_order)

## more parameters to be played with for different results -- try it!
nucl_min = 50
nucl_max = 150
thresh_max = 6

"""
Primary function of program:

This section of code loops through the normalized scores.
If a score is above 0, that position becomes the potential start position of a nucleosome.
If a following score drops below 0, the threshold counter is increased by 1. If there are 6 or more
consecutive scores below 0, the end point is recorded.

The start/end points are kept if the length between start/end is greater than nucl_min but less
than nucl_max base pairs long, bounds that specify the potential size of a nucleosome. This parameters can
be adjusted to grant varrying results in what is accepted as a potential nucleosome.

The threshold can also be adjusted (thresh_max).
"""
for line in smooth_list:
    if line > 0 and not started:
        started = True
        start_base = pos + startpos
          
    if line < 0 and started:
        thresh += 1
  
    if thresh >= thresh_max:
        thresh = 0
        started = False
        end_base = pos + startpos
        region_len = end_base - start_base
        if region_len >= nucl_min and region_len <= nucl_max:
            nucl.add_interval(Interval(start_base, end_base))
            count += 1

    pos += 1
## initialize a list for data from find_contig() method
contig_list = []

for read in nucl.find(0, end_base + 1):
    ## for each potential nucleosome interval, call find contig method.
    find_contig(read, smooth_list, startpos, contig_list)

## lastly, call this method to print all probable nucleosome ranges to a bed file
calc_nucl_distance(contig_list, smooth_list, startpos, chrom_num)
