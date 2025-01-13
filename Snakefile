import csv
import os
from datetime import datetime
import configparser
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

print(config)
modules = config['modules']

# ----- Pipeline Parameters -----
base_sample_sheet = f'user_input/{config["base_sample_sheet"]}'
countess_sample_ini = config["countess_sample_ini"]
barcode_variant_map = config["barcode_variant_map"]
count_cutoff = config["count_cutoff"]
sample_filter = config['sample_filter'] if bool(config['sample_filter']) else 'Complete'
bcl2fastq_argumets = config['bcl2fastq_arguments']

raw_data_directory = config['input_data_path']
# ----- Optional CutAdapt Parameter -----
cutadapt_trim = bool(config['cutadapt_trim_boolean'])

# ----- Run Configuration -----
run_path = config['run_name'] if config['run_name'] else config['run_timestamp']

threads = config['threads']
memory = config['memory']
# -----------------

if config['run_type'] == 'demux':
    step_output = [f'{run_path}/demux.txt', f'{run_path}/prep_fastqs.txt']
else:
    step_output = [f'{run_path}/countess.txt', f'{run_path}/metrics/correlation_matrix.jpeg']


rule all:
    input: step_output
    run:
        print('#----- Complete -----#')

# Only needed in some instances, place around the I7 and I5 sequences (index and index2) in samplesheet
# def index_reverse_compliment_and_trim(barcode_string, trim_amount):
#     reverse_compliment = {
#         'A': 'T',
#         'T': 'A',
#         'C': 'G',
#         'G': 'C'
#     }
#     reverse_compliment_string = ''
#     for letter in barcode_string:
#         reverse_compliment_string = reverse_compliment[letter] + reverse_compliment_string
    
#     reverse_compliment_string[:trim_amount]
#     return reverse_compliment_string[:-trim_amount]
        
rule clean_and_filter_samplesheet:
    input: base_sample_sheet
    output: f'{run_path}/{sample_filter}_samplesheet.csv'
    run:
        print('#----- Create Sample Sheet With Expected bcl2fastq Column Names -----#')

        bcl2fastq_samplesheet_rows = ['Sample_ID', 'Sample_Name', 'I7_Index_ID', 'index', 'I5_Index_ID', 'index2', 'Sample_Project']
        cleaned_samples = [['[Data]','','','','','',''], bcl2fastq_samplesheet_rows]

        optional_cutadapt_data = [['sample_id', 'f_amp_primer', 'r_amp_primer', 'f_amp_min_overlap', 'r_amp_min_overlap', 'target_length', 'error']] if cutadapt_trim == True else ''

        with open(f'{input}', newline='') as samplesheet:
            samples = csv.DictReader(samplesheet)
            
            for sample in samples:
                if sample_filter == 'Complete' or sample_filter in sample['Sample_ID'][:len(sample_filter)]:
                    clean_sample = [
                        sample['Sample_ID'].replace('.', '_'),
                        sample['Sample_Name'].replace('.', '_'),
                        sample[config['samplesheet_params']['i7_index_id_column_name']],
                        sample[config['samplesheet_params']['i7_index_sequence_column_name']],
                        sample[config['samplesheet_params']['i7_index_id_column_name']],
                        sample[config['samplesheet_params']['i5_index_sequence_column_name']],
                        sample['Sample_Project']
                    ]
                    cleaned_samples.append(clean_sample)
                    
                    if cutadapt_trim == True:
                        sample_to_trim = [
                            sample['Sample_ID'].replace('.', '_'),
                            sample[config['cutadapt_params']['f_amp_primer_column_name']].upper(),
                            sample[config['cutadapt_params']['r_amp_primer_column_name']].upper(),
                            config['cutadapt_params']['f_amp_min_overlap'],
                            config['cutadapt_params']['r_amp_min_overlap'],
                            sample[config['cutadapt_params']['target_length_column_name']],
                            config['cutadapt_params']['error']
                        ]
                        optional_cutadapt_data.append(sample_to_trim)

        if cutadapt_trim == True:
            print('----- Creating Trim Data File for CutAdapt -----')
            with open(f'{run_path}/{sample_filter}_trim_data.csv', 'w+', newline='', encoding='utf-8') as cutadapt_csv:
                csv.writer(cutadapt_csv).writerows(optional_cutadapt_data)

        with open(f'{output}', 'w+', newline='', encoding='utf-8') as cleaned_samplesheet:
            csv.writer(cleaned_samplesheet).writerows(cleaned_samples)


rule demux_and_pair:
    input: f'{run_path}/{sample_filter}_samplesheet.csv'
    output: f'{run_path}/demux.txt'
    params:
        input_data_dir = raw_data_directory,
        run_path = run_path,
        bcl2fastq_arguments = config['bcl2fastq_arguments'],
        processing_threads = threads - 4,
        bcl2fastq_loading_threads = 2,
        bcl2fastq_writing_threads = 2
    shell:
        """
        module load bcl2fastq/2.20
        module load pear/0.9.11

        pwd={workflow.basedir}
        samplesheet=$pwd/{input}

        cd $pwd/{params.run_path}
        
        mkdir -p ./bcl2fastq_output/
        echo {params.input_data_dir}

        bcl2fastq -R "{params.input_data_dir}" -o ./bcl2fastq_output/ --sample-sheet "$samplesheet" --no-lane-splitting -p {params.processing_threads} -r {params.bcl2fastq_loading_threads} -w {params.bcl2fastq_writing_threads} {params.bcl2fastq_arguments}

        mkdir -p ./pear_output/

        ls -d ./bcl2fastq_output/*/ | tr '\\n' '\\0' | xargs -0 -n 1 basename | grep -v -E 'Reports|Stats' | while read FOLDER; do
            echo "$FOLDER"
            mkdir ./pear_output/${{FOLDER}}
        
            ls -1 ./bcl2fastq_output/${{FOLDER}}/*_R1_001.fastq.gz | tr '\\n' '\\0' | xargs -0 -n 1 basename | sed 's/_R1_001.fastq.gz//' | while read SAMPLE; do
                echo "$SAMPLE"
                mkdir ./pear_output/${{FOLDER}}/${{SAMPLE}}
                pear -f ./bcl2fastq_output/${{FOLDER}}/${{SAMPLE}}_R1_001.fastq.gz -r ./bcl2fastq_output/${{FOLDER}}/${{SAMPLE}}_R2_001.fastq.gz -o ./pear_output/${{FOLDER}}/${{SAMPLE}}/${{SAMPLE}} -q 30 -j {params.processing_threads} -n 10 
        
            done
        
        done
        
        (echo "Complete") > {output}
        """
                              

# rule trim_fastqs:
#     input: f'{run_path}/{sample_filter}_samplesheet.csv'
#     run:
#         if cutadapt_trim == True:
#             print('#----- RUNNING CUTADAPT & FASTQC -----#')
            
#             pear_output_paths = get_fastq_file_paths(f'{input}', run_path)
    
#             samples_to_trim = []
    
#             with open(f'{run_path}/{sample_filter}_trim_data.csv') as trim_csv:
#                 samples_to_trim = list(csv.DictReader(trim_csv))
            
#             for pear_output_path in pear_output_paths:
#                 assembled_fastq_file = pear_output_path.split('/')[-1]
#                 assembled_fastq_path = f'{pear_output_path}/{assembled_fastq_file}'
                
#                 sample_dir_name = fastq_file.split('.')[0]
#                 sample_id = '_'.join(sample_dir_name.split('_')[:-1])
                
#                 trimmed_sample_data = [sample for sample in samples_to_trim if sample['sample_id'] == sample_id][0]
#                 f_amp_primer = trimmed_sample_data['f_amp_primer']
#                 r_amp_primer = trimmed_sample_data['r_amp_primer']
#                 f_amp_min_overlap = trimmed_sample_data['f_amp_min_overlap']
#                 r_amp_min_overlap = trimmed_sample_data['r_amp_min_overlap']
#                 target_length = trimmed_sample_data['target_length']
#                 error = trimmed_sample_data['error']
#                 trimmed_file_name = f'{sample_dir_name}.assembled.trimmed.fastq'
    
#                 shell(f'cutadapt -g \"{f_amp_primer};min_overlap={f_amp_min_overlap}...{r_amp_primer};min_overlap={r_amp_min_overlap}\" -e {error} -m {target_length} -M {target_length} -o {pear_output_path}/{trimmed_file_name} {pear_output_path}/{assembled_fastq_file}')

# def get_extensions(cutadapt_trim):
#     if cutadapt_trim == True:
#         return ['assembled.fastq', 'assembled.trimmed.fastq', 'discarded.fastq', 'unassembled.forward.fastq', 'unassembled.reverse.fastq']
#     else:
#         return ['assembled.fastq', 'discarded.fastq', 'unassembled.forward.fastq', 'unassembled.reverse.fastq']


def get_fastq_file_paths(samplesheet, run_path):
    try:
        pear_output_path = f'{workflow.basedir}/{run_path}/pear_output'
        all_fastq_directories = []
    
        list_of_sample_projects = os.listdir(pear_output_path)
        for sample_project in list_of_sample_projects:
            all_fastq_directories = all_fastq_directories + os.listdir(f'{pear_output_path}/{sample_project}')
            
        list_of_sample_paths = []
        with open(samplesheet, newline='') as samplesheet:
            next(samplesheet)
            samples = csv.DictReader(samplesheet)
            for sample in samples:
                sample_dir_name = [sample_name for sample_name in all_fastq_directories if sample['Sample_ID'] in sample_name]
                if len(sample_dir_name) > 0:
                    sample_project = sample['Sample_Project']
                    list_of_sample_paths.append(f'{sample_project}/{sample_dir_name[0]}/{sample_dir_name[0]}')
        print(f'list paths{list_of_sample_paths}')
        return list_of_sample_paths 
    except:
        return []
                       
# ----- Run CutAdapt and FastQC -----
rule prep_fastqs_for_countess:
    input: f'{run_path}/demux.txt'
    output: f'{run_path}/prep_fastqs.txt'
    params:
        run_path = run_path,
        fastq_paths = get_fastq_file_paths(f'{run_path}/{sample_filter}_samplesheet.csv', run_path),
        fastqc_module = modules['fastqc'],
    shell:
        """
        FASTQC_OUTPUT_PATH="{params.run_path}/fastqc_output"
        mkdir -p $FASTQC_OUTPUT_PATH

        module load {params.fastqc_module}

        for FASTQ_PATH in {params.fastq_paths}; do
            FOLDER=$(basename $(dirname $FASTQ_PATH))
            FASTQC_THREADS=$(find $FASTQ_PATH/. -name '*.fastq' | wc -l)
            
            mkdir -p $FASTQC_OUTPUT_PATH/$FOLDER
            mkdir -p $FASTQC_OUTPUT_PATH/$FOLDER/$(basename $FASTQ_PATH)

            echo "fastqc $FASTQ_PATH/*.fastq -t $FASTQC_THREADS --outdir $FASTQC_OUTPUT_PATH/$FOLDER/$(basename $FASTQ_PATH)"
            fastqc $FASTQ_PATH/*.fastq -t $FASTQC_THREADS --outdir $FASTQC_OUTPUT_PATH/$FOLDER/$(basename $FASTQ_PATH)
            
        done
        (echo "Complete") > {output}
        """

# ----- Run CountESS -----
rule run_countess_vampseq:
    input: f'user_input/{countess_sample_ini}'
    output: f'{run_path}/countess.txt'
    run:
        print('#----- RUNNING COUNTESS -----#')
        
        # countess_ini = configparser.ConfigParser()
        # countess_ini.read(f'{input}')

        # countess_ini['LoadBarcodeMap']['files.0.filename'] = f"'{barcode_variant_map}'"
        # countess_ini['SaveTotalCounts']['filename'] = f"'../{run_path}/unfiltered_counts.csv'"
        # countess_ini['SaveCountBeforeFreq']['filename'] = f"'../{run_path}/programmed_count.csv'"
        # countess_ini['SaveTotalFreq']['filename'] = f"'../{run_path}/programmed_freq.csv'"
        # countess_ini['SaveScores']['filename'] = f"'../{run_path}/final_scores.csv'"
        # countess_ini['FilterLowCountVars']['filters.0.columns.0.value'] = f"'{count_cutoff}'"

        
        # with open(f'{input}', 'w') as modified_ini:    # save
        #     countess_ini.write(modified_ini)
                    
        shell(f'countess_cmd {input} && (echo "Complete") > {output}')

# ----- Create Plots -----
rule create_plots:
    input: f'{run_path}/countess.txt'
    output: f'{run_path}/metrics/correlation_matrix.jpeg'
    run:
        print('#----- Creating Plots -----#')
        shell(f'mkdir -p {run_path}/metrics')
        unfiltered_counts = pd.read_csv(f'{run_path}/unfiltered_counts.csv')
        num_count_cols = unfiltered_counts.columns.tolist()[1:]
        total_counts = unfiltered_counts[num_count_cols].sum()
        unique_barcodes = unfiltered_counts[num_count_cols].astype(bool).sum(axis = 0)
        shell(f'mkdir -p {run_path}/metrics/read_hists')

        def plot_histogram(count_series, sample):
            # plt.figure(figsize=(8, 6))
            plt.hist(count_series, histtype='step', log=True, bins=10**np.linspace(0, 4.5, 100))
            plt.xscale('log')
            plt.yscale('log')
            plt.title(f'Read Count Histogram - {sample}')
            plt.xlabel('Read Counts (Log Scale)')
            plt.ylabel('Frequency (Log Scale)')
            plt.grid(True)
            plt.savefig(f'{run_path}/metrics/read_hists/{sample}.jpeg')
            
        for col in num_count_cols:
            plot_histogram(unfiltered_counts[col], col)

        vampseq_output = pd.read_csv(f'../{run_path}/final_scores.csv')
        fig, axs = plt.subplots(1, 3)
        rep = 1
        for ax in axs:
            ax.hist(vampseq_output[vampseq_output['set'] == 'Missense'][f'score_rep_{rep}'], bins = 30, color = 'lightgrey', histtype = 'bar', label = 'Missense', edgecolor = 'black')
            ax.hist(vampseq_output[vampseq_output['set'] == 'Synonymous'][f'score_rep_{rep}'], bins = 30, color = 'blue', histtype = 'bar', label = 'Synonymous', edgecolor = 'black', alpha = 0.7)
            ax.hist(vampseq_output[vampseq_output['set'] == 'Nonsense'][f'score_rep_{rep}'], bins = 30, color = 'red', histtype = 'bar', label = 'Nonsense', edgecolor = 'black', alpha = 0.5)
            ax.legend()
            ax.set_xlabel('VAMPseq scores')
            ax.set_ylabel('# of Variants')
            fig.set_figwidth(20)
            rep += 1
        fig.savefig(f'{run_path}/metrics/stacked_score_count_hist.jpeg')

                       
        counts_and_freq = pd.read_csv(f'{run_path}/programmed_freq.csv')
        for i in range(1, 4):
            counts_and_freq[f'sum_rep_{i}'] = counts_and_freq[f'count__bin_1__rep_{i}'] + counts_and_freq[f'count__bin_2__rep_{i}'] + counts_and_freq[f'count__bin_3__rep_{i}'] + counts_and_freq[f'count__bin_4__rep_{i}']
            
        rep_counts = counts_and_freq[['aaChanges'] + list(counts_and_freq.columns[-3:])]
        joined = vampseq_output.set_index('aaChanges').join(rep_counts.set_index('aaChanges'))

        fig, axs = plt.subplots(1, 3)
        rep = 1
        for ax in axs:
            x = joined[f'sum_rep_{rep}']
            y = joined[f'score_rep_{rep}']
            ax.scatter(x, y, s = 5, alpha = 0.5, edgecolor = 'black', linewidth = 0.15)
            ax.set_xscale('log')
            ax.set_xlabel('Count (Log Scale)')
            ax.set_ylabel('VAMPseq score')
            fig.set_figwidth(20)
            rep += 1
        fig.savefig(f'{run_path}/metrics/score_total_count_scatter.jpeg')

        
        columns_of_interest = ['score_rep_1', 'score_rep_2', 'score_rep_3']  # Replace with actual column names
        # Select these columns
        df_selected = vampseq_output[columns_of_interest]
        # Create a PairGrid object and map scatter plots to the upper triangle
        g = sns.PairGrid(df_selected)
        g = g.map_lower(sns.scatterplot)
        # Map histograms or density plots to the diagonal
        g = g.map_diag(sns.histplot, kde=True)
        # Compute the correlation matrix
        correlation_matrix = df_selected.corr()
        # Display the correlation matrix
        # Adjust the plot
        plt.subplots_adjust(top=0.95)
        plt.savefig(f'{output}')
