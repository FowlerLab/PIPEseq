import subprocess
import csv
import os
import datetime
import configparser

# ----- Pipeline Parameters -----
base_sample_sheet = config['base_sample_sheet']
countess_sample_ini = config['countess_sample_ini']
raw_data_directory = config['raw_data_directory']
sample_filter = config['sample_filter'] if bool(config['sample_filter']) else 'Complete'

# ----- Optional CutAdapt Parameter -----
cutadapt_trim = bool(config['cutadapt_trim'] == True)

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


        
rule clean_and_filter_samplesheet:
    input: base_sample_sheet
    output: f'{run_path}/{sample_filter}_samplesheet.csv'
    run:
        print('#----- Create Sample Sheet With Expected bcl2fastq Column Names -----#')

        bcl2fastq_samplesheet_rows = ['Sample_ID', 'Sample_Name', 'I7_Index_ID', 'I5_Index_ID', 'index', 'index2', 'Sample_Project']
        cleaned_samples = [['[Data]','','','','','',''], bcl2fastq_samplesheet_rows]

        optional_cutadapt_data = [['sample_id', 'f_amp_primer', 'r_amp_primer', 'f_amp_min_overlap', 'r_amp_min_overlap', 'target_length', 'error']] if cutadapt_trim == True else ''

        with open(f'{input}', newline='') as samplesheet:
            samples = csv.DictReader(samplesheet)
            
            for sample in samples:
                if sample_filter == 'Complete' or sample_filter in sample['Sample_ID'][:len(sample_filter)]:
                    clean_sample = [
                        sample['Sample_ID'].replace('.', '_'),
                        sample['Sample_Name'].replace('.', '_'),
                        sample['I7_Index_ID'],
                        sample['I5_Index_ID'],
                        sample['I7_Index_seq'],
                        sample['I5_Index_seq'],
                        sample['Sample_Project']
                    ]
                    cleaned_samples.append(clean_sample)
                    
                    if cutadapt_trim == True:
                        sample_to_trim = [
                            sample['Sample_ID'].replace('.', '_'),
                            sample[config['cutadapt_f_amp_primer_column']].upper(),
                            sample[config['cutadapt_r_amp_primer_column']].upper(),
                            config['cutadapt_f_amp_min_overlap'],
                            config['cutadapt_r_amp_min_overlap'],
                            sample[config['cutadapt_target_length_column']],
                            config['cutadapt_error']
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
    params:
        input_data_dir = raw_data_directory,
        run_path = run_path,
        bcl2fastq_processing_threads = threads - 4,
        bcl2fastq_loading_threads = 2,
        bcl2fastq_writing_threads = 2
    shell:
        """
        module load bcl2fastq/2.20
        module load pear/0.9.11

        pwd={workflow.basedir}
        raw_data_directory=$pwd/raw_data/{params.input_data_dir}
        samplesheet=$pwd/{input}

        cd $pwd/{params.run_path}
        
        mkdir ./bcl2fastq_output/
                
        bcl2fastq -R $raw_data_directory -o ./bcl2fastq_output/ --sample-sheet $samplesheet --no-lane-splitting -p {params.bcl2fastq_processing_threads} -r {params.bcl2fastq_loading_threads} -w {params.bcl2fastq_writing_threads}
                
        mkdir ./pear_output/

        ls -d ./bcl2fastq_output/*/ | tr '\\n' '\\0' | xargs -0 -n 1 basename | grep -v -E 'Reports|Stats' | while read FOLDER; do
            echo "$FOLDER"
            mkdir ./pear_output/${{FOLDER}}
        
            ls -1 ./bcl2fastq_output/${{FOLDER}}/*_R1_001.fastq.gz | tr '\\n' '\\0' | xargs -0 -n 1 basename | sed 's/_R1_001.fastq.gz//' | while read SAMPLE; do
                echo "$SAMPLE"
                mkdir ./pear_output/${{FOLDER}}/${{SAMPLE}}
                pear -f ./bcl2fastq_output/${{FOLDER}}/${{SAMPLE}}_R1_001.fastq.gz -r ./bcl2fastq_output/${{FOLDER}}/${{SAMPLE}}_R2_001.fastq.gz -o ./pear_output/${{FOLDER}}/${{SAMPLE}}/${{SAMPLE}} -q 3 -j {params.bcl2fastq_processing_threads} -n 10 --barcode-mismatch 0
        
            done
        
        done
        """


# ----- Run CutAdapt and FastQC -----
rule prep_fastqs_for_vampseq:
    params:
        samplesheet = 'run-240807_VH00979_267_AAFTYGJM5-202408190629-Complete/Complete_samplesheet.csv',
        run_path = 'run-240807_VH00979_267_AAFTYGJM5-202408190629-Complete',
        cutadapt_trim = cutadapt_trim
    run:
        print('#----- RUNNING CUTADAPT & FASTQC -----#')
        shell(f'mkdir -p {params.run_path}/fastqc_output')
        
        pear_output_paths = get_fastq_file_paths('run-240807_VH00979_267_AAFTYGJM5-202408190629-Complete/Complete_samplesheet.csv', params.run_path)

        list_of_fastqs_for_vampseq = []
        samples_to_trim = []

        if params.cutadapt_trim == True:
            with open(f'{params.run_path}/{sample_filter}_trim_data.csv') as trim_csv:
                samples_to_trim = list(csv.DictReader(trim_csv))

        for pear_output_path in pear_output_paths:
            output_fastqs = os.listdir(f'{workflow.basedir}/{pear_output_path}')
            for fastq_file in output_fastqs:
                shell(f'fastqc {pear_output_path}/{fastq_file} --outdir {params.run_path}/fastqc_output')
                if '.assembled' in fastq_file:
                    if params.cutadapt_trim == True:
                        sample_dir_name = fastq_file.split('.')[0]
                        sample_id = '_'.join(sample_dir_name.split('_')[:-1])
                        
                        trimmed_sample_data = [sample for sample in samples_to_trim if sample['sample_id'] == sample_id][0]
                        f_amp_primer = trimmed_sample_data['f_amp_primer']
                        r_amp_primer = trimmed_sample_data['r_amp_primer']
                        f_amp_min_overlap = trimmed_sample_data['f_amp_min_overlap']
                        r_amp_min_overlap = trimmed_sample_data['r_amp_min_overlap']
                        target_length = trimmed_sample_data['target_length']
                        error = trimmed_sample_data['error']
                        trimmed_file_name = f'{sample_dir_name}.assembled.trimmed.fastq'
            
                        shell(f'cutadapt -g \"{f_amp_primer};min_overlap={f_amp_min_overlap}...{r_amp_primer};min_overlap={r_amp_min_overlap}\" -e {error} -m {target_length} -M {target_length} -o {pear_output_path}/{trimmed_file_name} {pear_output_path}/{fastq_file}')
                        
                        shell(f'fastqc {workflow.basedir}/{pear_output_path}/{trimmed_file_name} --outdir {workflow.basedir}/{params.run_path}/fastqc_output')
                        list_of_fastqs_for_vampseq.append(f'{pear_output_path}/{trimmed_file_name}')
                    
                    else:
                        list_of_fastqs_for_vampseq.append(f'{pear_output_path}/{fastq_file}')
                        
        print(f'run{list_of_fastqs_for_vampseq}')
                

----- Run CountESS -----v
rule run_countess_vampseq:
    params:
        run_path = 'run-240807_VH00979_267_AAFTYGJM5-202408190629-Complete',
        cutadapt_trim = cutadapt_trim
    run:
        print('#----- RUNNING COUNTESS -----#')
        shell(f'mkdir {params.run_path}/countess_inis')

        pear_output_paths = get_fastq_file_paths('run-240807_VH00979_267_AAFTYGJM5-202408190629-Complete/Complete_samplesheet.csv', params.run_path)

        files_for_countess = []

        for pear_output_path in pear_output_paths:
            sample_name = pear_output_path.split('/')[-1]
            if params.cutadapt_trim == True:
                files_for_countess.append(f'{pear_output_path}/{sample_name}.assembled.trimmed.fastq')
            else:
                files_for_countess.append(f'{pear_output_path}/{sample_name}.assembled.fastq')
                
        print(f'count{files_for_countess}')
    
        
        default_vampseq_ini = 'vampseq_default.ini'
        default_vampseq_config = configparser.ConfigParser()
        default_vampseq_config.read(default_vampseq_ini)
        
        print(default_vampseq_config.sections())
        
        ini_file_name = f'{params.run_path}/countess_inis/{run_timestamp}_{sample_filter}.ini'

        auto_countess_config = configparser.ConfigParser()

        for key in default_vampseq_config:
            if key != 'DEFAULT':
                print(key)
                auto_countess_config[key] = {}
                for sub_key in default_vampseq_config[key]:
                    if key == 'FASTQ Load' and 'filename' in sub_key:
                        for index in range(len(files_for_countess)):
                            auto_countess_config['FASTQ Load'][f'files.{index}.filename'] = f"'../{files_for_countess[index]}'"
                    if key == 'Save Scores' and 'filename' in sub_key:
                            auto_countess_config['Save Scores']['filename'] = f"'../{sample_filter}_vampseq.csv'"
                    else:
                        auto_countess_config[key][sub_key] = default_vampseq_config[key][sub_key]
        with open(ini_file_name, 'w') as vampseq_ini:
            auto_countess_config.write(vampseq_ini)

        # shell(f'countess_cmd {ini_file_name}')
