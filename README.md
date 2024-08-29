# PIPEseq

### A modular pipeline for VAMPseq and LABALseq (soon&trade;)

#### Pre-Configuration:
- conda (https://conda.io/projects/conda/en/latest/user-guide/install/index.html)
- bcl2fastq/2.20 (in cluster modules, technically bcl2fastq2)
- pear/0.9.11 (in cluster modules)
- fastqc/0.12.1 (in cluster modules)

#### Python Dependencies (will be installed from yaml or Setup)
- Python 3.10
- CutAdapt 4.9 (bioconda)
- CountESS (pip)
  
#### Setup
- Create a Conda environment with yaml file `conda env create -f conda_env.yaml -n <env_name>`
- Activate the Conda environment `conda activate <env_name>`
- Install CountESS with pip `pip install countess` (depending on Python/Pip configuration, may use pip3)



### The Pipeline

