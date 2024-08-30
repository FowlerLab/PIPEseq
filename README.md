# PIPEseq

### A modular pipeline for VAMPseq and LABALseq (soon&trade;)

#### Pre-Configuration:
- conda (https://conda.io/projects/conda/en/latest/user-guide/install/index.html)
- bcl2fastq/2.20 (in cluster modules, technically bcl2fastq2)
- pear/0.9.11 (in cluster modules)
- fastqc/0.12.1 (in cluster modules)

#### Python Dependencies (will be installed from yaml or Setup)
- Python 3.10
- CutAdapt 4.9
- Snakemake 7.32.4
- CountESS (via pip)
  
#### Setup
- Clone PIPEseq into your directory
- In your terminal, navigate to the PIPEseq directory, `cd <directory_with_sequence_data>/PIPEseq`
- Create a Conda environment with yaml file `conda env create -f conda_env.yaml -n <env_name>`
- Activate the Conda environment `conda activate <env_name>`
- Copy your Samplesheet for bcl2fastq and barcode-variant map for VAMP-seq into the PIPEseq directory

#### Expected Folder Structure
- *<directory_with_raw_data>/
    - PIPEseq/
        - conda_env.yaml
        - README.md
        - run.sh
        - Snakefile
        - *<samplesheet_file>.csv
        - *<countess_file>.ini

*=your file (or folder)


### The Pipeline
- Rules:
    - clean_and_filter_samplesheet
    - demux_and_pair
    - prep_fastqs_for_countess
    - 