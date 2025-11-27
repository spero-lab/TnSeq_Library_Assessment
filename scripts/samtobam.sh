#!/bin/bash

# Activate conda env ( w/ yq & fastqc) 
    conda activate TNSeq_env

# MANAGING VARIABLES / INPUTS / OUTPUTS 
Config_Template=../config/config_template.yml
# Parse .yml --> isolate whether mapq filter is being applied or not 
MAPQ_filter=$(yq eval '.adjustable_parameters.MAPQ_Filter.MAPQ_filter' "$Config_Template")
MAPQ_threshold=$(yq eval '.adjustable_parameters.MAPQ_Filter.MAPQ_threshold' "$Config_Template")
sample_name=$(yq eval '.sample_name' "$Config_Template") #used to help ID newly cat-ed .fastq.gz 

# 1. Convert SAM --> sorted BAM 
samtools view -bS ../intermediate/"${sample_name}.sam" | samtools sort -o ../intermediate/"${sample_name}.sorted.bam"

# 2. MAPQ Filter (op.) + #3 Calc. % Mapped Reads 
    # If No, cont. (skip)
    if [[ "$MAPQ_filter" == "no" ]]; then
        echo "Skipping mapq filtering"
        # Calc. % mapped reads --> append to "${sample_name}_performance_summary.tsv 
            # Append metric name (for .tsv) 
        printf "%s\t" "% mapped reads (no mapq filter)" >> ../intermediate/"${sample_name}_performance_summary.tsv"  #%s\t is just saying to add a \t after the proceding string
            # Append value (for .tsv) 
        # COMMENTS FOR PIPE BELOW 
        #flagstat outputs a summary on perf. stats of the bam file (including % mapped reads)to terminal
            # isolating only the line (from flagstat) that lists the % mapped reads
            # a. Breaks line to be delimited by (,%, and ). b. Grabbing the 2nd field ($2) isolates the # in b/w "(" and "%". [ ex. line --> 1000 + 0 mapped ( 85.00% : N/A)] — Note: awk automatically adds a \n at end of isolated #. 
            # appends this as an entry to the "${sample_name}_performance_summary.tsv file 
        samtools flagstat ../intermediate/"${sample_name}.sorted.bam" \
            | grep "mapped (" \
            | awk -F '[]()%[]' '{print $2}' >> ../intermediate/"${sample_name}_performance_summary.tsv"
    ## If Yes, then apply mapq filter 
    elif [[ "$MAPQ_filter" == "yes" ]]; then
        echo "applying mapq filter" 
        # (A) Calc. % mapped reads pre-filter --> append to "${sample_name}_performance_summary.tsv 
            # Append metric name (for .tsv) 
        printf "%s\t" "% mapped reads (pre-mapq-filter)" >> ../intermediate/"${sample_name}_performance_summary.tsv" 
            # Append value (for .tsv) 
        samtools flagstat ../intermediate/"${sample_name}.sorted.bam" \
            | grep "mapped (" \
            | awk -F '[]()%[]' '{print $2}' \
            >> ../intermediate/"${sample_name}_performance_summary.tsv"
        # (B) Apply MAPQ Filter 
        samtools view -b -q $MAPQ_threshold ../intermediate/"${sample_name}.sorted.bam" > ../intermediate/"${sample_name}.filtered.bam"
        # (C) Calc. % mapped reads post-filter --> append to "${sample_name}_performance_summary.tsv 
        # Append metric name (for .tsv) 
        printf "%s\t" "% mapped reads (post-mapq-filter)" >> ../intermediate/"${sample_name}_performance_summary.tsv"  #%s\t is just saying to add a \t after the proceding string 
            # Append value (for .tsv) 
        samtools flagstat ../intermediate/"${sample_name}.filtered.bam" \
            | grep "mapped (" \
            | awk -F '[]()%[]' '{print $2}' \
            >> ../intermediate/"${sample_name}_performance_summary.tsv" 
    else
        echo "Invalid input. Please answer yes or no in config_template.yml's MAPQ_filter parameter" 
    fi 
exit
