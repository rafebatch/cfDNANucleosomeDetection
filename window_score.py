import time
import os
import sys
import math
import pysam
from bx.intervals.intersection import Intersecter, Interval
from optparse import OptionParser, OptionGroup

"""
TO RUN: 
$ python window_score.py BAM_file.bam -f output_file.wig

You can also specifiy a region to look at with -r:
$ python window_score.py BAM_file.bam -r $REGION -f output_file.wig

Just make sure $REGION is in the format -- chrom:start-end (e.g. 1:250000-300000)
Make sure your bam_file is sorted/indexed. Not sure if the pipe option works.
"""

## prints out the time for running at the end -- your welcome
start_time = time.time()

"""
Option handling
"""
parser = OptionParser("%prog [options] BAMfile")
parser.add_option("-p", "--PIPE", dest="pipe", help="Read BAMfile into program from PIPE, write back to PIPE", default = False, action="store_true")
parser.add_option("-o", "--outdir", dest="outdir", help="Create output files in another directory.", default = None)
parser.add_option("-r","--region", dest="region", help="Region to be looked up. If no region specified, region will be defined by start/end reads of file passed in.")
parser.add_option("-f", "--file", dest="filename", help="File to be created storing output data.", default=None)

(options, args) = parser.parse_args()   #load option values in options, leftover args in args

"""
Handle input file
"""
if options.pipe:            ## if file piped from samtools, specify option -p
    infile = pysam.Samfile('-', "rb")
elif len(args) == 0:       ## if no file supplied, error raised
    sys.stderr.write("Please supply a BAM file to read from or specify PIPE option\n")
    sys.exit()
else:                       ## if not piped, read file from command line
    filename = args[0]
    infile = pysam.Samfile(filename, "rb")

"""
Handle output files/pipe
"""
if options.filename != None:
    outfilename = options.filename
else:
    sys.stderr.write("No output file name specified - saving to default_output_score.wig\n")
    outfilename = "default_output_score.wig"

if options.pipe:                    ## if -p specified, pass file along pipe without writing
    outfile = sys.stdout
elif options.outdir != None:        ## otherwise, write to output file
    outfilepath = os.path.join(os.getcwd(), options.outdir + "/" + outfilename)
    directory = os.path.dirname(outfilepath)        ## get the directory full path
    
    if not os.path.exists(directory):               ## check if it exists
        os.makedirs(directory)                      ## if not, create the directory
    outfile = open(outfilepath, "w+")               ## write the file to the directory
else:
    outfile = open(outfilename, "w+")               ## otherwise, write to current dir

"""
Handle data
"""
total_reads = 0
start_end_list = Intersecter()      ## list to hold start/end of each read in infile

init_pos = 0        ## for defining region if user does not
final_pos = 0       ## for defining region if user does not

read_start = 0      ## start of each read
read_end = 0        ## end of each read
other_read_end = 0


"""
Main functionality:
"""

for read in infile:
    total_reads += 1
    
    ## get start position of the first read
    if total_reads == 1:
        init_pos = read.reference_start
    
    ## if read is not a proper pair or is a fail, duplicate, etc, continue
    if read.is_qcfail or read.is_duplicate or read.qual == None:
        continue
    
    ## only want proper pair reads
    if read.is_proper_pair:         ## if read is proper pair, add entire PE read start/end to list
        if abs(read.template_length) < 10000:       ## precaution for weird reads - sometimes have huge length
            if read.is_read1:
                
                ## get start/end of read
                read_start = min(read.reference_start, read.next_reference_start) + 1
                read_end = read_start + abs(read.template_length)
                
                ## add start/end to interval list
                start_end_list.add_interval(Interval(read_start, read_end))

## get last position from last read
final_pos = min(read.reference_start, read.next_reference_start) + abs(read.template_length)

if options.region != None:
    chrom = options.region.split(':')[0]
    start, end = map(int, options.region.split(':')[1].split('-'))
else:
    chrom = read.reference_name
    start, end = init_pos, final_pos

window_size = 120               ## define the window size (could be made an option)
prot_region = window_size//2    ## definitely a parameter worth messing with

## write the first line of info for the peak calling file
outfile.write("fixedStep chrom=chr%s start=%d step=1" % (chrom, start))

## This is where the scoring occurs
for pos in range(start,end+1):
    ## get the start/end points of the window in current position
    w_start, w_end = pos - prot_region, pos + prot_region  ## define the window at each position in region
    end_points = 0                                         ## set to 0 each window shift
    intact_reads = 0
   
    ## get the number of fragments within the window that are fully intact (+1)
    ## or have an endpoint (-1)
    for read in start_end_list.find(w_start, w_end):
        if (read.start > w_start) or (read.end < w_end):
            end_points += 1
        else:
            intact_reads += 1
    ## write the score for each window position to the output wig file
    outfile.write("\n" + str(intact_reads - end_points))        

## close the file
outfile.close()

## print total execution time
sys.stderr.write("--- %s seconds ---\n" % (time.time() - start_time))
