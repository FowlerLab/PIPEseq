#!/bin/bash

# ----- Parameters (Optional) -----
RAW_DATA_DIRECTORY=
THREADS=24
MEMORY=36
# ----- End Parameters -----

RUN_TIMESTAMP=$(date +%Y%m%d_%H%M) 

# snakemake \
#     --cores 4 \
#     --latency-wait 120 \
#     --configfile "user_variables.yaml" \
#     --config threads=$THREADS memory=$MEMORY \
#     run_timestamp=$RUN_TIMESTAMP raw_data_directory=${RAW_DATA_DIRECTORY:-$(basename $(cd ../ && pwd))}
    
snakemake \
    --cluster "qsub -l mfree=${MEMORY}G -pe serial $THREADS" \
    --cluster-cancel "qdel" \
    --jobs 10 \
    --latency-wait 120 \
    --configfile "user_variables.yaml" \
    --config threads=$THREADS memory=$MEMORY \
    run_timestamp=$RUN_TIMESTAMP raw_data_directory=${RAW_DATA_DIRECTORY:-$(basename $(cd ../ && pwd))}
