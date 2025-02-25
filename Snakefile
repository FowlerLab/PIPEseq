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
barcode_variant_map = f'user_input/{config["barcode_variant_map"]}'
count_ini = f'default_files/{config["count_ini"]}'
score_ini = f'user_input/{config["score_ini"]}'
default_cutoff = config["default_cutoff"]
sample_filter = config['sample_filter'] if bool(config['sample_filter']) else 'Complete'
bcl2fastq_argumets = config['bcl2fastq_arguments']
raw_data_directory = config['input_data_path']

# ----- Run Configuration -----
run_path = config['run_name'] if config['run_name'] else config['run_timestamp']
threads = config['threads']
memory = config['memory']

# ----- Optional CutAdapt Parameter -----
cutadapt_trim = bool(config['cutadapt_trim_boolean'])
# -----------------

if config['run_type'] == 'demux':
    step_output = f'{run_path}/prep_fastqs.txt'
elif config['run_type'] == 'count':
    step_output = [f'{run_path}/prep_fastqs.txt', f'{run_path}/metrics/raw_count_hists.jpeg']
else:
    step_output = [f'{run_path}/metrics/correlation_matrix.jpeg']
    # step_output = [f'{run_path}/countess_score.csv', f'{run_path}/metrics/correlation_matrix.jpeg']


rule all:
    input: step_output
    run:
        print('#----- Complete -----#')


# ----- Convert samplesheet to the one bcl2fastq expects
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


# -----  Demultiplex and Pair Illumina sequences -----
rule demux_and_pair:
    input: f'{run_path}/{sample_filter}_samplesheet.csv'
    output: temp(f'{run_path}/demux.txt')
    params:
        input_data_dir = raw_data_directory,
        run_path = run_path,
        processing_threads = config['threads'],
        bcl2fastq_module = modules['bcl2fastq'],
        pear_module = modules['pear'],
        bcl2fastq_arguments = config['bcl2fastq_arguments']
    shell:
        """
        module load {params.bcl2fastq_module}
        module load {params.pear_module}
        pwd={workflow.basedir}
        samplesheet=$pwd/{input}
        (echo "Complete") > {output}
        cd $pwd/{params.run_path}
        
        mkdir -p ./bcl2fastq_output/
        echo {params.input_data_dir}

        bcl2fastq -R "{params.input_data_dir}" -o ./bcl2fastq_output/ --sample-sheet "$samplesheet" --no-lane-splitting {params.bcl2fastq_arguments}

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
        """


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
        fastqc_module = modules['fastqc']
    shell:
        """
        FASTQC_OUTPUT_PATH="{params.run_path}/fastqc_output"
        mkdir -p $FASTQC_OUTPUT_PATH

        module load {params.fastqc_module}

        for FASTQ_PATH in {params.fastq_paths}; do
            FOLDER=$(basename $(dirname $FASTQ_PATH))
            FASTQC_THREADS=$(find {params.run_path}/pear_output/$FASTQ_PATH* -name '*.fastq' | wc -l)
            echo "$FASTQC_OUTPUT_PATH/$FOLDER"
            mkdir -p $FASTQC_OUTPUT_PATH/$FOLDER

            echo "fastqc {params.run_path}/pear_output/$FASTQ_PATH*.fastq -t $FASTQC_THREADS --outdir $FASTQC_OUTPUT_PATH/$FOLDER/$(basename $FASTQ_PATH)"
            fastqc {params.run_path}/pear_output/$FASTQ_PATH*.fastq -t $FASTQC_THREADS --outdir $FASTQC_OUTPUT_PATH/$(basename $FASTQ_PATH)
            
        done
        (echo "Complete") > {output}
        """


#----- Get raw counts -----
rule run_countess_counting:
    input: expand(f'{run_path}/pear_output/{{fastq}}.assembled.fastq', fastq = get_fastq_file_paths(f'{run_path}/{sample_filter}_samplesheet.csv', run_path))
    output: f'{run_path}/cutoffs.csv'
    conda: 'rule_envs/countess_env.yaml'
    params:
        run_path = run_path,
        count_ini = count_ini,
        default_cutoff = default_cutoff
    shell:
        '''
        mkdir -p {params.run_path}/counts
        (echo "filename,count") > {output}
        
        for FASTQ_PATH in {input}; do
            FILENAME=$(basename $FASTQ_PATH)
            SAMPLE_ID=${{FILENAME%.assembled*}}
            (echo "$SAMPLE_ID,{params.default_cutoff}") >> {output}
            echo "$FASTQ_PATH"
            echo countess_cmd --set LoadFASTQ.files.0.filename='"'$FASTQ_PATH'"' --set SaveRawCounts.filename='"'{params.run_path}/counts/$SAMPLE_ID.csv'"' {params.count_ini}
            countess_cmd --set LoadFASTQ.files.0.filename='"'$FASTQ_PATH'"' --set SaveRawCounts.filename='"'{params.run_path}/counts/$SAMPLE_ID.csv'"' --log=debug {params.count_ini}
        done
        '''

                       
# ----- Create cutoff histograms -----
rule create_cutoff_hists:
    input: f'{run_path}/cutoffs.csv'
    output: f'{run_path}/metrics/raw_count_hists.jpeg'
    params: 
        samples = [f'{run_path}/counts/{x.split('/')[-1]}.csv' for x in get_fastq_file_paths(f'{run_path}/{sample_filter}_samplesheet.csv', run_path)]
    run:
        shell(f'mkdir -p {run_path}/metrics')
        count_files = params.samples

        reps = int(len(count_files) / 4)
        fig, axes = plt.subplots(reps, 4, figsize=(25,15))
        
        total_count_obj = {}
        unique_barcode_obj = {}
        
        x = 0
        y = 0
        count_files.sort()
        for count_file in count_files:
            sample_id = count_file.split('/')[-1][:-4]
            unfiltered_counts = pd.read_csv(f'{count_file}')
            total_count_obj[sample_id] = unfiltered_counts['count'].sum()
            unique_barcode_obj[sample_id] = len(unfiltered_counts)
            
            ax = axes[y, x]
            ax.hist(unfiltered_counts, histtype='step', log=True, bins=10**np.linspace(0, 4.5, 100))
            ax.set_title(f'Read Count Histogram - {sample_id}')
            ax.set_xlabel('Read Counts (Log Scale)')
            ax.set_ylabel('Frequency (Log Scale)')
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.grid(True)
            
            if x == 3:
                x = 0
                y += 1
            else:
                x += 1
        
        fig.tight_layout() 
        fig.savefig(f'{output}')

        # Read Count
        fig, axes = plt.subplots(2, 1, figsize=(15,15))
        axes[0].bar(list(total_count_obj.keys()), list(total_count_obj.values()), color='g')
        axes[0].tick_params(axis = 'x', labelrotation = 90)
        axes[0].set_title('Total Reads')
        
        #unique read count
        axes[1].bar(list(unique_barcode_obj.keys()), list(unique_barcode_obj.values()), color='b')
        axes[1].tick_params(axis = 'x', labelrotation = 90)
        axes[1].set_title('Unique Reads')
        fig.tight_layout()
        fig.savefig(f'{run_path}/metrics/total_and_unique_reads.jpeg')

           
# ----- Run Scoring INI file -----
rule run_countess_scoring:
    input: f'{run_path}/cutoffs.csv'
    output: f'{run_path}/countess_score.csv'
    conda: 'rule_envs/countess_env.yaml'
    params:
        run_path = run_path,
        score_ini = score_ini,
        barcode_variant_map = barcode_variant_map,
    shell:
        '''
        echo countess_cmd --set LoadCutoff.files.0.filename='"'{input}'"' --set LoadBarcodeMap.files.0.filename='"'{params.barcode_variant_map}'"' --set "'"LoadTotalCount.files.0.filename='"'{params.run_path}/counts/*.csv'"'"'" --set SaveScores.filename='"'{output}'"' --log=debug {params.score_ini}
        countess_cmd --set LoadCutoff.files.0.filename='"'{input}'"' --set LoadBarcodeMap.files.0.filename='"'{params.barcode_variant_map}'"' --set "'"LoadTotalCount.files.0.filename='"'{params.run_path}/counts/*.csv'"'"'" --set SaveScores.filename='"'{output}'"' --log=debug {params.score_ini}
        '''

           
# ----- Create Plots -----
# rule create_plots:
#     input: f'{run_path}/countess_score.csv'
#     output: f'{run_path}/metrics/correlation_matrix.jpeg'
#     run:
#         print('#----- Creating Plots -----#')

#         countess_output = pd.read_csv(f'{input}')
rule create_plots:
    output: f'{run_path}/metrics/correlation_matrix.jpeg'
    run:
        print('#----- Creating Plots -----#')

        countess_output = pd.read_csv(f'{run_path}/countess_score.csv')
        fig, axes = plt.subplots(1, 3, figsize=(15,5))
        
        # measurements = ['score__rep_1', 'score__rep_2', 'score__rep_3']
        measurements = ['score_rep_1', 'score_rep_2', 'score_rep_3']
        
        for ax, meas in zip(axes, measurements):
            sns.kdeplot(data=countess_output, x=meas, hue='set', fill=True, ax=ax)
            ax.set_title(meas)
        
        # plt.tight_layout(rect=[0, 0, 1, 0.95])  # To ensure the title fits nicely
        fig.savefig(f'{run_path}/metrics/stacked_hist.jpeg')

                       
        for i in range(1, 3):
            countess_output[f'sum_rep_{i}'] = countess_output[f'count__bin_1__rep_{i}'] + countess_output[f'count__bin_2__rep_{i}'] + countess_output[f'count__bin_3__rep_{i}']

        fig, axs = plt.subplots(1, 3)
        rep = 1
        for ax in axs:
            x = countess_output[f'sum_rep_{rep}']
            y = countess_output[f'score_rep_{rep}']
            ax.scatter(x, y, s = 5, alpha = 0.5, edgecolor = 'black', linewidth = 0.15)
            ax.set_xscale('log')
            ax.set_xlabel('Count (Log Scale)')
            ax.set_ylabel('VAMPseq score')
            fig.set_figwidth(20)
            rep += 1
        fig.savefig(f'{run_path}/metrics/score_total_count_scatter.jpeg')

        # Correlation Matrix
        columns_of_interest = ['score_rep_1', 'score_rep_2', 'score_rep_3']  # Replace with actual column names
        # Select these columns
        df_selected = countess_output[columns_of_interest]
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


# -------------------- CutAdapt/one-off code to clean --------------------
                          
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
              

# ----- CutAdapt Trimming -----
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