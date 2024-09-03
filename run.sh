#!/bin/bash

# ----- Parameters (Optional) -----
RAW_DATA_DIRECTORY=
THREADS=12
MEMORY=12
# ----- End Parameters -----

# module load fastqc/0.12.1

# snakemake \
#     --cores 1 \
#     --latency-wait 120 \
#     --configfile "user_variables.yaml" \
#     --config threads=$THREADS memory=$MEMORY \
#     run_timestamp=$(date +%Y%m%d_%H%M) raw_data_directory=${RAW_DATA_DIRECTORY:-$(basename $(cd ../ && pwd))}
    
snakemake \
    --cluster "qsub -l mfree=${MEMORY}G -pe serial $THREADS" \
    --cluster-cancel "qdel" \
    --jobs 10 \
    --latency-wait 120 \
    --configfile "user_variables.yaml" \
    --config threads=$THREADS memory=$MEMORY \
    run_timestamp=$(date +%Y%m%d_%H%M) raw_data_directory=${RAW_DATA_DIRECTORY:-$(basename $(cd ../ && pwd))}
