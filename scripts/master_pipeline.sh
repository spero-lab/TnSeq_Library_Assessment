!/bin/bash

# Changing file permissions for all .sh scripts in scripts folder (in case not already done) 
chmod 755 ./* 

# Running Scripts 
    # Set up Conda Env  
    ./setup1_conda.sh
    # Cp Input Files 
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate TNSeq_env
    ./setup2_organize.sh
    # Run FastQC (pre-filter)
    ./fastqc_pre.sh
    # Filter Reads 
        # Activate conda env ( w/ yq) 
            # conda activate TNSeq_env — already done in previous line 
        # Run Filter Reads Script — will need to use YQ to grab the input file paths / parameters from config_template.yml to input as arguments when running the script (argparse)
            # Manage Variables / Input / Outputs 
            Config_Template=../config/config_template.yml
            # Parse .yml --> isolate custom QC parameters + filepath names 
            tbs=$(yq eval '.adjustable_parameters.QC.transposon_border_sequence' "$Config_Template")
            hdt=$(yq eval '.adjustable_parameters.QC.hamming_distance_threshold_tbs' "$Config_Template")
            ps=$(yq eval '.adjustable_parameters.QC.phredscore_threshold' "$Config_Template")
            hdq=$(yq eval '.adjustable_parameters.QC.hamming_distance_threshold_qs' "$Config_Template")
            sample_name=$(yq eval '.sample_name' "$Config_Template") #used to help ID newly cat-ed .fastq.gz 
            # Run filter_reads.py script 
            ./filter_reads.py -f ../data/"${sample_name}.fastq.gz" \
                -o ../intermediate/"${sample_name}_filtered_cleaned.fastq" \
                -m ../intermediate/"${sample_name}_filtered.fastq" \
                -u ../intermediate/"${sample_name}_removed.fastq" \
                -s ../intermediate/"${sample_name}_performance_summary.tsv" \
                -t "${tbs}" \
                -ht "${hdt}" \
                -q "${ps}" \
                -hq "${hdq}" \
    # Run FastQC (post-filter)
    ./fastqc_post.sh
    # Index Reference Genome (.fna) — Increases Alignment Speed 
        # Parse .yml --> isolate reference genome strain name 
        ref_organism=$(yq eval '.reference_genome_organism' "$Config_Template")
        # Run Bowtie2 
        bowtie2-build -f ../data/*.fna "$ref_organism"
        # Move output .bt2 files from scripts/ --> intermediate/  
        mv *.bt2 ../intermediate
    # Alignment (End-to-End) 
        # Run Alignment Script via Bowtie2 
        ./alignment.sh
    # SAM --> BAM 
        ./samtobam.sh
    
    # Pipeline is still under development. 


        # Deactivate conda env (w/ yq)
        # conda deactivate 
        # X (op.) — Zip the output fastq files to save space? 
exit