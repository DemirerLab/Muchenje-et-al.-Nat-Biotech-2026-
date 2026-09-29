# R2 25S integration site mapping
The workflow outlined here is used for analyzing ONT amplicon long-read sequencing data of R2 retrotransposon insertion sites in plants to map each insertion to one or multiple 25S sites within the genome.

## Set up
This workflow has been developed and tested on MacOS Apple Silcon. Of note is the use of `barrnap` to identify 25S sites, this requires HMMER 3.1b2 which functions on x64 architecture. Therefore for the barrnap process we set the architecture to x64 and use Rosetta. Please install Rosetta if not already installed `softwareupdate --install-rosetta`. 
We provide two environment files, _r2_environment.yaml_ for everything except the barrnap process, and _barrnap_environment.yaml_ for running barrnap.

### R2 Environment Install
```bash
conda env create -f r2_environment.yaml
```

### barrnap Environment Install
```bash
CONDA_SUBDIR=osx-64 conda env create -f barrnap_environment.yaml
```

Please use the R2 environment for all steps except `barrnap 25S annotation`. 

## Step 1: File organization
For the workflow you require a genome fasta file, in our analysis we use the _Nicotiana benthamiana_ genome from https://doi.org/10.1038/s41477-024-01849-y. This file can be found at [CNCB](https://ngdc.cncb.ac.cn/) under project number PRJCA022857. You also require at least one fastq file from the ONT sequencing, our four replicate _N. benthamiana_ fastq files are available at [SRA](https://www.ncbi.nlm.nih.gov/sra) under project number **PRJNA1463879**. To be compatible with the workflow your directory should look like this:
```
- R2-25S-integration-mapping-main (name of this folder can be changed)
│   ├── workflow.ipynb
│   ├── data
│   │   └── ONT_A.fastq
│   │   └── ONT_B.fastq
│   │   └── ONT_C.fastq
│   │   └── ONT_D.fastq
│   ├── reference
│   │   └── genome.fasta
```
## Step 2: barrnap 25S annotation
In your R2-25S-integration-mapping-main directory, first activate the barrnap environment and then run barrnap in eukaryotic mode. The output is a genome.gff file with the locations of each predicted rDNA sequence from your genome file. For the 2.88 Gb _N. benthamiana_ genome this takes ~30 minutes with 10 threads on M4 Pro. Adjust thread count as desired.
```bash
conda activate barrnap_environment
conda config --env --set subdir osx-64
barrnap --kingdom euk --threads 10 reference/genome.fasta  > reference/genome.gff
```
## Step 3: Running the analysis pipeline
Switch conda environments for all other steps.
```bash
conda activate r2_environment
```
Open the workflow.ipynb Jupyter Notebook and follow the in-built instructions. All output files and plots matching those in our publication will be generated and saved.
