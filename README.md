**RNA Editing Analysis Toolkit**  
This repository contains scripts for analyzing A-to-I RNA editing from both short-read (NGS) and long-read (Nanopore) sequencing data.


**NGS_cut_bam.py：**  
Python script for cutting reads in BAM files with multiple length options and direction controls  
Usage:python NGS_cut_bam.py [-h] -i INPUT --prefix PREFIX --length LENGTH  
  [--direct {left,random,right}] [--paired]


**NGS_cut_FASTQ.py：**  
Simple FASTQ read cutter for quality trimming or length adjustment  
Usage:python NGS_cut_FASTQ.py input.fastq output.fastq length or python NGS_cut_FASTQ.py input.fastq.gz output.fastq length


**LRS_cut_bam.py：**  
Enhanced BAM cutter with additional functionality for Nanopore data  
Usage:python LRS_cut_bam.py [-h] -i INPUT --prefix PREFIX --length LENGTH  
  [--direct {left,random,right}] [--paired]


**LRS_cut_FASTQ_window.py：**  
Python script for sliding window-based FASTQ cutting  
Usage:python LRS_cut_FASTQ_window.py input.fastq output.fastq window_size step_size
