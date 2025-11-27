#!/bin/bash

# Set up Conda Env. / Install Necessary Packages 

# TNSeq_env 
    # Miniconda Installation: User should have miniconda installed beforehand 
    # Install all packages specified in .yaml 
        # Source / locate user's conda 
        source "$(conda info --base)/etc/profile.d/conda.sh"
        # create & activate conda env 
        conda env create -f ../config/environment.yml
        conda activate TNSeq_env
    # UNIT-TEST 
        echo "TNSeq_env pkg versions"
        # Check package versions 
        python --version
        fastqc --version
        bowtie2 --version
        samtools --version
        yq --version

    # switching conda envs to setup Multiqc_env 
    conda deactivate 
# Multiqc_env
    # create & activate conda env 
        conda create --name Multiqc_env python=3.8
        conda activate Multiqc_env
    # Install MultiQC 
        #conda install -c bioconda multiqc=1.30 # conda couldn't solve this package dependency 
        pip install multiqc==1.30
    # UNIT-TEST 
        echo "MultiQC_env pkg versions"
        # Check package versions 
        python --version
        multiqc --version
exit
