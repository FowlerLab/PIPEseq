#!/bin/bash

# ----- Parameters (Optional) -----
THREADS=32
MEMORY=64

# ----- End Parameters -----
    
snakemake \
    --executor cluster-generic \
    --cluster-generic-submit-cmd "qsub -l mfree=${MEMORY}G -l h=fl003" \
    --cluster-generic-cancel-cmd "qdel" \
    --software-deployment-method conda \
    --jobs 30 \
    --latency-wait 120 \
    --use-envmodules \
    --configfile "user_variables.yaml" \
    --config run_timestamp=$(date +%Y%m%d_%H%M) \
    threads=$THREADS memory=$MEMORY
