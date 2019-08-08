# Nucleosome Detection and Mapping using cell-free DNA from Blood Plasma

This repository contains code for the detection and position mapping of nucleosomes from sequenced blood plasma.

## Window Scores
The "window_score.py" file takes a BAM file as input and outputs a WIG file containing scores for each base pair position along the reference genome, for the chromosome specified by the BAM file.

Currently, window_score.py only supports one chromosome at a time, so use samtools, i.e.
```bash
samtools view -b in.bam 1 > in_chr1.bam
```
would generate a new bam file with only data from chromosome 1.

window_score.py can then be called:
```bash
python window_score.py input_BAM_chr1.bam -f output_chr1.wig
```

Output is stored in output_chr1.wig. Additional options can be found within window_score.py

## Peak Calling

The peak_call.py method then takes a WIG file produced by window_score.py and outputs a mapping consisting of the probable locations of nucleosomes in BED format.

To run peak_call.py:
```bash
python peak_call.py input_wig.wig > output_bed.bed
```

Please ensure that you specify a wig file as the input.

## Feature Extract

This program simply outputs the number of nucleosomes in a user defined bin size. This can be used for analysis in the comparison of nucleosomes per region between samples.

Usage:

```bash
python feature_extract.py example_chr2.bed -s unhealthy -c breast -w 1000000 > example_chr2_ana.txt
```

This would result in a text file containing the number of nucleosomes per 1000000 bp bin. The -s tag is required to indicate a sample is unhealthy, the -c tag is to indicate which cancer type the sample is afflicted by. If healthy, specify -s healthy and -c none.

## Viewing Results

Both the BED and WIG files produced from the original BAM file can be viewed with the IGV program. This program nicely plots the nucleosomes, as well as the peaks in the WIG file. Can be useful for visualization of data.

 
