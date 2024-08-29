#!/bin/bash

# ----- Pipeline Parameters -----
BASE_SAMPLE_SHEET=SignalPeptides_sampleSheet_3May24run.csv
COUNTESS_SAMPLE_INI=vampseq_default.ini
RAW_DATA_DIRECTORY=
SAMPLE_FILTER=

# ----- Cutadapt Parameters (Optional) -----
CUTADAPT_TRIM_BOOLEAN=True
CUTADAPT_ERROR=0.1
CUTADAPT_F_AMP_MIN_OVERLAP=3
CUTADAPT_R_AMP_MIN_OVERLAP=3

CUTADAPT_F_AMP_PRIMER_COLUMN_NAME=upstream_flanking_seq
CUTADAPT_R_AMP_PRIMER_COLUMN_NAME=downstream_flanking_seq
CUTADAPT_TARGET_LENGTH_COLUMN_NAME=Target_length

# ----- Cluster Parameters (Optional) -----
THREADS=12
MEMORY=12

# ----- End Parameters -----

snakemake -R demux_and_pair \
    --cluster "qsub -l mfree=${MEMORY}G -pe serial $THREADS" \
    --cluster-cancel "qdel" \
    --jobs 10 \
    --latency-wait 120 \
    --config base_sample_sheet=${BASE_SAMPLE_SHEET:-$(ls -1 | grep 'ample')} countess_sample_ini=${COUNTESS_SAMPLE_INI:-$(ls -1 | grep '.ini')} \
    raw_data_directory=${RAW_DATA_DIRECTORY:-$(ls raw_data)} sample_filter=${SAMPLE_FILTER:-False} \
    cutadapt_trim=${CUTADAPT_TRIM_BOOLEAN:-False} cutadapt_error=${CUTADAPT_ERROR:-False} \
    cutadapt_f_amp_min_overlap=${CUTADAPT_F_AMP_MIN_OVERLAP:-False} cutadapt_r_amp_min_overlap=${CUTADAPT_R_AMP_MIN_OVERLAP:-False} \
    cutadapt_f_amp_primer_column=${CUTADAPT_F_AMP_PRIMER_COLUMN_NAME:-False} cutadapt_r_amp_primer_column=${CUTADAPT_R_AMP_PRIMER_COLUMN_NAME:-False} \
    cutadapt_target_length_column=${CUTADAPT_TARGET_LENGTH_COLUMN_NAME:-False} \
    run_timestamp=$(date +%Y%m%d_%H%M) threads=$THREADS memory=$MEMORY

# module load fastqc/0.12.1

# snakemake -R run_countess_vampseq \
#     --cores 1 \
#     --config base_sample_sheet=${BASE_SAMPLE_SHEET:-$(ls -1 | grep 'ample')} countess_sample_ini=${COUNTESS_SAMPLE_INI:-$(ls -1 | grep '.ini')} \
#     raw_data_directory=${RAW_DATA_DIRECTORY:-$(ls raw_data)} sample_filter=${SAMPLE_FILTER:-False} \
#     cutadapt_trim=${CUTADAPT_TRIM_BOOLEAN:-False} cutadapt_error=${CUTADAPT_ERROR:-False} \
#     cutadapt_f_amp_min_overlap=${CUTADAPT_F_AMP_MIN_OVERLAP:-False} cutadapt_r_amp_min_overlap=${CUTADAPT_R_AMP_MIN_OVERLAP:-False} \
#     cutadapt_f_amp_primer_column=${CUTADAPT_F_AMP_PRIMER_COLUMN_NAME:-False} cutadapt_r_amp_primer_column=${CUTADAPT_R_AMP_PRIMER_COLUMN_NAME:-False} \
#     cutadapt_target_length_column=${CUTADAPT_TARGET_LENGTH_COLUMN_NAME:-False} \
#     run_timestamp=$(date +%Y%m%d_%H%M) threads=$THREADS memory=$MEMORY