**RNA Editing Analysis Toolkit**  
This repository contains scripts for analyzing A-to-I RNA editing from both short-read (NGS) and long-read (Nanopore) sequencing data.


**NGS_cut_bam.py：**  
Python script for cutting reads in BAM files with multiple length options and direction controls.  
Usage: python NGS_cut_bam.py [-h] -i INPUT --prefix PREFIX --length LENGTH  
  [--direct {left,random,right}] [--paired]


**NGS_cut_FASTQ.py：**  
Simple FASTQ read cutter for quality trimming or length adjustment.  
Usage: python NGS_cut_FASTQ.py input.fastq output.fastq length or python NGS_cut_FASTQ.py input.fastq.gz output.fastq length


**LRS_cut_bam.py：**  
Enhanced BAM cutter with additional functionality for Nanopore data.  
Usage: python LRS_cut_bam.py [-h] -i INPUT --prefix PREFIX --length LENGTH  
  [--direct {left,random,right}] [--paired]


**LRS_cut_FASTQ_window.py：**  
Python script for sliding window-based FASTQ cutting.  
Usage: python LRS_cut_FASTQ_window.py input.fastq output.fastq window_size step_size

RNA_Edit_Analysis
Index files: BAM files must be indexed (.bam.bai or .bai) before running any script. The wrapper does not create indices automatically.
Memory usage: For large site lists (e.g., 3 million sites), samtools mpileup may use significant memory. We recommend limiting to ≤5 million sites per run.
Base quality offset: The scripts assume Phred+33 encoding (standard for modern FASTQ). If your data uses Phred+64, change $offset in Query_Editing_Level.pl to 64.
Parallelism: For many samples, use the provided parallel example or a job scheduler (e.g., SGE, SLURM).

1. Batch editing level calculation (NGS data, Illumina)
Use the wrapper editing_site_call.pl to process all BAM files in a directory.
The script expects indexed BAMs (.bam.bai or .bai).
Default minimum base quality = 30 (suitable for Illumina).
Usage:perl editing_site_call.pl --in /path/to/bam_directory --out /path/to/output_directory
--in : directory containing .bam files (and their index files)
--out : output directory where results will be stored

2. Single‑sample editing level calculation (Nanopore / LRS data)
For long‑read data (e.g., Oxford Nanopore), use Query_Editing_Level.pl directly.
Specify the minimum base quality as the 4th argument (default = 30, set to 7 for Nanopore).
Usage:perl Query_Editing_Level.pl <sites.bed> <sample.bam> <output.txt> [min_base_qual]
Example (Nanopore, quality ≥7):perl Query_Editing_Level.pl neg_control_3M.bed sample.bam sample.all.levels.txt 7
Example (Illumina, quality ≥30 – can omit the last argument):perl Query_Editing_Level.pl known_sites.bed sample.bam sample.all.levels.txt

3. Editing level classification (Alu / non‑Alu / other)
After obtaining .all.levels.txt files for all samples, run the classification script:
Usage:bash level_analysis_adapt.sh /path/to/call_level_output_directory
Input directory should contain subdirectories (one per sample) each with a .all.levels.txt file.
The script uses bedtools intersect to assign each site to one of three categories:
ALU : overlaps an Alu repeat
nonALU_rep : overlaps a non‑Alu repeat
other : does not overlap any repeat

Example Workflow
1.Prepare negative control sites (e.g., 3 million random non‑edited A sites).
2.Run editing level calculation for each dataset:
NGS: perl editing_site_call.pl --in ./bams/ --out ./results/
Nanopore: perl Query_Editing_Level.pl sites.bed sample.bam sample.out.txt 7
3.Classify results:bash level_analysis_adapt.sh ./results/
4.Analyze the classified files (e.g., plot editing levels per category).
