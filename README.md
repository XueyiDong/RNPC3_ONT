# Using long- and short-read RNA-seq to study minor class intron splicing

This repository contains code used to perform the analysis and generate the figures for my PhD thesis chapter 4.

## Files and folders

[aligned_minimap2](aligned_minimap2): genomic alignment BAM files and scripts for mapping

[bambu](bambu): scripts and outputs for *bambu*

<!-- [count](count): scripts and outputs from running *featureCounts* for gene-level counting -->

<!-- [data](data): a soft link to the original data -->

<!-- [figures](figures): figures organizaed to be used in thesis -->

[salmon](salmon): scripts and outputs for *salmon*

[short_bam_merged](short_bam_merged): merge bams of samples from the same group for visualization of coverage

[transcript_analysis](transcript_analysis): scripts and outputs for transcriptomic analysis: quality control, DTE, DTU, visualization of read coverage

## DTE and DTU results table

These can be found in [transcript_analysis/results](transcript_analysis/results) (for long read) and [transcript_analysis/results_short](transcript_analysis/results_short) (for short read). 