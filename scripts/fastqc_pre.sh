#!/bin/bash

# Activate conda env ( w/ yq & fastqc) 
    conda activate TNSeq_env

# MANAGING VARIABLES / INPUTS / OUTPUTS 
Config_Template=../config/config_template.yml
    # Parse .yml --> PREFIX in name of cat-ed FASTQ file
    Sample_Name=$(yq '.sample_name' "$Config_Template")

# Pre-Filter FASTQC 
    # Create folder to store fastq results 
    mkdir ../results/fastqc
    # Run FastQC on the concatenated fastq file 
    fastqc ../data/"$Sample_Name".fastq.gz -o ../results/fastqc
