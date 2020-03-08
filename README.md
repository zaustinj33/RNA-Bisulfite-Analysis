# RNA-Bisulfite-Analysis
An in-depth exploration of RNA Bisulfite data in R

## 1) Filter by RNAseq data, and create replicate overlap and union files
Script for generating filtered overlap and union files for meRanCall output files
Input: 2 biological replicate files labelled *_[0-9]_m5c* stored in <home>/RNA_BS/m5C_siteLists/
        2 RNAseq gene counts files with the same beginning text as the meRanCall files
        optional 2 RNAseq transcript files 
        RNAseq files stored in <home>/RNAseq/gene_lists

 Output: 2 RNAseq-filtered m5C counts files
          One overlap file
          One union file
          All files contain: coverage, count, methylation level, and seq context info

The user only needs to specify home_dir and project. As long as files are in the
right locations, and files and columns are labelled correctly this will work.
There are no checkpoints to protect the user from these problems, so please hastag 
and write functions until all tests have passed.

## 2) Visualization of replicate overlap data
### a) Dot plot and Venn diagram of count data
We first plot comparisons between replicates to ensure consistency between samples
### b) Methylation level histogram and sequence context
We then look at average methylation level and sequence context for initial impressions of the data

If data is deemed 'good' then we move further downstream to 2) site location analysis and 3) differential methylation analysis

## 3) Site location preparation and analysis


