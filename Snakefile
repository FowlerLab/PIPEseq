import csv
import os
from datetime import datetime
import configparser

print(config)

# ----- Pipeline Parameters -----
base_sample_sheet = f'user_input/{config["base_sample_sheet"]}'
countess_sample_ini = f'user_input/{config["countess_sample_ini"]}'
barcode_variant_map = f'user_input/{config["barcode_variant_map"]}'
sample_filter = config['sample_filter'] if bool(config['sample_filter']) else 'Complete'
bcl2fastq_argumets = config['bcl2fastq_arguments']

raw_data_directory = config['raw_data_directory']
# ----- Optional CutAdapt Parameter -----
cutadapt_trim = bool(config['cutadapt_trim_boolean'])

# ----- Run Configuration -----
run_timestamp = config['run_timestamp']
run_name = f'{run_timestamp}-{sample_filter}'
run_path = f'run-{raw_data_directory}-{run_name}'

threads = config['threads']
demux_and_pair_threads = int((threads-1) * .75)
memory = config['memory']
# -----------------


def get_fastq_file_paths(samplesheet, run_path):
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
                list_of_sample_paths.append(f'{run_path}/pear_output/{sample_project}/{sample_dir_name[0]}')
    print(f'list paths{list_of_sample_paths}')
    return list_of_sample_paths


if config['steps'] == 1:
    step_output = f'{run_path}/demux.txt'
else:
    step_output = f'{run_path}/countess_inis/{run_timestamp}_{sample_filter}.ini'


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
        bcl2fastq_processing_threads = threads - 4,
        bcl2fastq_loading_threads = 2,
        bcl2fastq_writing_threads = 2
    shell:
        """
        module load bcl2fastq/2.20
        module load pear/0.9.11

        pwd={workflow.basedir}
        raw_data_directory=$pwd/../../{params.input_data_dir}
        samplesheet=$pwd/{input}

        (echo "Complete") >  {output}
        cd $pwd/{params.run_path}
        
        mkdir ./bcl2fastq_output/
        echo $raw_data_directory

        echo "bcl2fastq -R $raw_data_directory -o ./bcl2fastq_output/ --sample-sheet $samplesheet --no-lane-splitting -p {params.bcl2fastq_processing_threads} -r {params.bcl2fastq_loading_threads} -w {params.bcl2fastq_writing_threads} {params.bcl2fastq_arguments}"
        bcl2fastq -R $raw_data_directory -o ./bcl2fastq_output/ --sample-sheet $samplesheet --no-lane-splitting -p {params.bcl2fastq_processing_threads} -r {params.bcl2fastq_loading_threads} -w {params.bcl2fastq_writing_threads} {params.bcl2fastq_arguments}
                
        mkdir ./pear_output/

        ls -d ./bcl2fastq_output/*/ | tr '\\n' '\\0' | xargs -0 -n 1 basename | grep -v -E 'Reports|Stats' | while read FOLDER; do
            echo "$FOLDER"
            mkdir ./pear_output/${{FOLDER}}
        
            ls -1 ./bcl2fastq_output/${{FOLDER}}/*_R1_001.fastq.gz | tr '\\n' '\\0' | xargs -0 -n 1 basename | sed 's/_R1_001.fastq.gz//' | while read SAMPLE; do
                echo "$SAMPLE"
                mkdir ./pear_output/${{FOLDER}}/${{SAMPLE}}
                pear -f ./bcl2fastq_output/${{FOLDER}}/${{SAMPLE}}_R1_001.fastq.gz -r ./bcl2fastq_output/${{FOLDER}}/${{SAMPLE}}_R2_001.fastq.gz -o ./pear_output/${{FOLDER}}/${{SAMPLE}}/${{SAMPLE}} -q 30 -j {params.bcl2fastq_processing_threads} -n 10
        
            done
        
        done
        """



# ----- Run CutAdapt and FastQC -----
# rule prep_fastqs_for_countess:
#     input: f'{run_path}/{sample_filter}_samplesheet.csv'
#     params:
#         cutadapt_trim = cutadapt_trim
#     run:
#         print('#----- RUNNING CUTADAPT & FASTQC -----#')
#         shell(f'mkdir -p {run_path}/fastqc_output')
        
#         pear_output_paths = get_fastq_file_paths(input, run_path)

#         list_of_fastqs_for_vampseq = []
#         samples_to_trim = []

#         if params.cutadapt_trim == True:
#             with open(f'{params.run_path}/{sample_filter}_trim_data.csv') as trim_csv:
#                 samples_to_trim = list(csv.DictReader(trim_csv))

#         for pear_output_path in pear_output_paths:
#             output_fastqs = os.listdir(f'{workflow.basedir}/{pear_output_path}')
#             for fastq_file in output_fastqs:
#                 # shell(f'fastqc {pear_output_path}/{fastq_file} --outdir {run_path}/fastqc_output')
#                 if '.assembled' in fastq_file:
#                     if params.cutadapt_trim == True:
#                         sample_dir_name = fastq_file.split('.')[0]
#                         sample_id = '_'.join(sample_dir_name.split('_')[:-1])
                        
#                         trimmed_sample_data = [sample for sample in samples_to_trim if sample['sample_id'] == sample_id][0]
#                         f_amp_primer = trimmed_sample_data['f_amp_primer']
#                         r_amp_primer = trimmed_sample_data['r_amp_primer']
#                         f_amp_min_overlap = trimmed_sample_data['f_amp_min_overlap']
#                         r_amp_min_overlap = trimmed_sample_data['r_amp_min_overlap']
#                         target_length = trimmed_sample_data['target_length']
#                         error = trimmed_sample_data['error']
#                         trimmed_file_name = f'{sample_dir_name}.assembled.trimmed.fastq'
            
#                         shell(f'cutadapt -g \"{f_amp_primer};min_overlap={f_amp_min_overlap}...{r_amp_primer};min_overlap={r_amp_min_overlap}\" -e {error} -m {target_length} -M {target_length} -o {pear_output_path}/{trimmed_file_name} {pear_output_path}/{fastq_file}')
                        
#                         shell(f'fastqc {workflow.basedir}/{pear_output_path}/{trimmed_file_name} --outdir {workflow.basedir}/{run_path}/fastqc_output')
#                         list_of_fastqs_for_vampseq.append(f'{pear_output_path}/{trimmed_file_name}')
                    
#                     else:
#                         list_of_fastqs_for_vampseq.append(f'{pear_output_path}/{fastq_file}')
                        
#         print(f'{list_of_fastqs_for_vampseq}')
                

# ----- Run CountESS -----
rule run_countess_vampseq:
    input: f'{run_path}/demux.txt'
    output: f'{run_path}/countess_inis/{run_timestamp}_{sample_filter}.ini'
    run:
        print('#----- RUNNING COUNTESS -----#')
        shell(f'mkdir -p {run_path}/countess_inis')

        pear_output_paths = get_fastq_file_paths(f'{run_path}/{sample_filter}_samplesheet.csv', run_path)

        files_for_countess = []

        for pear_output_path in pear_output_paths:
            sample_name = pear_output_path.split('/')[-1]
            if cutadapt_trim == True:
                files_for_countess.append(f'{pear_output_path}/{sample_name}.assembled.trimmed.fastq')
            else:
                files_for_countess.append(f'{pear_output_path}/{sample_name}.assembled.fastq')
                
        print(f'count{files_for_countess}')
    
        
        default_vampseq_config = configparser.ConfigParser()
        default_vampseq_config.read(countess_sample_ini)
        
        print(default_vampseq_config.sections())
        print('hiii')
    

        auto_countess_config = configparser.ConfigParser()

        for key in default_vampseq_config:
            if key != 'DEFAULT':
                print(key)
                auto_countess_config[key] = {}
                for sub_key in default_vampseq_config[key]:
                    if key == 'FASTQ Load' and 'filename' in sub_key:
                        for index in range(len(files_for_countess)):
                            auto_countess_config['FASTQ Load'][f'files.{index}.filename'] = f"'../../{files_for_countess[index]}'"
                    if key == 'Barcode Map Load' and 'filename' in sub_key:
                        auto_countess_config['Barcode Map Load']['files.0.filename'] = f"'../../{barcode_variant_map}'"
                    if key == 'Save Scores' and 'filename' in sub_key:
                        auto_countess_config['Save Scores']['filename'] = f"'../{sample_filter}_countess_output.csv'"
                    else:
                        auto_countess_config[key][sub_key] = default_vampseq_config[key][sub_key]
        print('ahhh')
        print(auto_countess_config.sections())
        with open(f'{output}', 'w') as countess_ini:
            auto_countess_config.write(countess_ini)
            
        shell(f'countess_cmd {output}')
