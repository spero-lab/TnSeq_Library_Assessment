#!/bin/bash

# Activate conda env ( w/ yq & fastqc) 
    conda activate TNSeq_env

# MANAGING VARIABLES / INPUTS / OUTPUTS 
Config_Template=../config/config_template.yml
sample_name=$(yq eval '.sample_name' "$Config_Template") #used to help ID newly cat-ed .fastq.gz 
ref_organism=$(yq eval '.reference_genome_organism' "$Config_Template")

# ALIGNMENT 
    # Run Bowtie2 (Generate SAM file in intermediate/)
bowtie2 -x ../intermediate/$ref_organism \
    -U ../intermediate/"${sample_name}_filtered_cleaned.fastq" \
    -S ../intermediate/"${sample_name}.sam"

