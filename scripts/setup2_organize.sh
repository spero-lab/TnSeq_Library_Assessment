#!/bin/bash

# setup1.sh Pipeline_Folder_Name

# MANAGING VARIABLES / INPUTS / OUTPUTS 
Config_Template=../config/config_template.yml

# Move Input Files (location specified in .yaml) to Data Folder - at this point, file organization should have been conserved from Github. 
    # Parse .yml --> isolate pathway to all input files 
    Input_Fastq_Origin=$(yq '.input_files.all.input_fastq' "$Config_Template")
    Ref_Genome_Origin=$(yq '.input_files.all.input_ref' "$Config_Template")
    # Change permissions of Input Files 
    chmod 755 "$Input_Fastq_Origin"/*
    chmod 755 "$Ref_Genome_Origin"/*
    # cp input files --> Pipeline's Data folder 
    cp -r "$Input_Fastq_Origin"/* ../data
    cp -r "$Ref_Genome_Origin"/* ../data

# Cat FASTQ Files (if not already)
    # Parse .yml --> isolate whether fastq files are cat-ed or not & sample name 
    Fastq_cat=$(yq '.adjustable_parameters.input_data.fastq_cat' "$Config_Template")
    Sample_Name=$(yq '.sample_name' "$Config_Template")
    # If Yes, cont. (skip)
    if [[ "$Fastq_cat" == "yes" ]]; then
        echo "Skipping fastq concatenation"
    ## If no, then cat. fastq together (essentially, each pipeline run is for only ONE SAMPLE at a time, no per CONDITION)
    elif [[ "$Fastq_cat" == "no" ]]; then
        echo "Concatenating FASTQ Files" 
        # Merge all fastqs (different lanes/chunks), WHILE PRESERVING CORRECT NUMERIC ORDER (_001,_002, etc.)
			# sort -V ensures that files are processed in numeric order (lane order = L001, L002). 
			# xargs avoids errors such as 'argument list too long' since these files are so large, actually allows cl to apply a action to files previously isolated in the pipe 
		ls ../data/*.fastq.gz | sort -V | xargs zcat | gzip -c > ../data/"$Sample_Name".fastq.gz
    else
        echo "Invalid input. Please answer yes or no in config_template.yml" 
    fi 
exit