#!/bin/bash

# Activate conda env ( w/ yq & fastqc) 
    conda activate TNSeq_env

# MANAGING VARIABLES / INPUTS / OUTPUTS 
Config_Template=../config/config_template.yml
    # Parse .yml --> PREFIX in name of cat-ed FASTQ file
    Sample_Name=$(yq eval '.sample_name' "$Config_Template")

# Post-Filter FASTQC 
    # Run FastQC on the filtered (X cleaned, so transposon seq.s should still be in this one) fastq file 
    fastqc ../intermediate/"$Sample_Name"_filtered.fastq -o ../results/fastqc