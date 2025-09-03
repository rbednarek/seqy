# Seq-Y
`/'siː.ki/`

## NGS Analysis Pipelines
Pipeline written in Nextflow and Snakemake for standard NGS analysis. Includes read processing, alignment and optional UMI-based deduplication functionality.  

## Nextflow

### Step 1: Clone the Repository

```bash
git clone https://github.com/rbednarek/seqy.git
cd seqy
```


### Step 2: Set Up Conda Environment

Run the config script to create a new conda environment and install all required dependencies.

```bash
bash ./config.sh nextflow
```

This script will:
- Create a new conda environment named `seqy`
- Install nextflow and other required software packages
- Note: If you want to install both snakemake and nextflow in this environment just run `bash ./config.sh`

### Step 3: Activate the Envrionment

After the setup script completes, activate the conda environment if necessary:

```bash
conda activate seqy
```

### Step 4: Run Pipeline

Basic Usage:

```bash
nextflow run ./nextflow_pipeline/seqy.nf --samplesheet samples.csv --reference_genome genome.fa
```
    
Required parameters:
--samplesheet       CSV file with columns: sample_name,R1,R2
--reference_genome  Path to reference genome FASTA file

    
Optional parameters:
--outdir           Output directory (default: ./results)
--merged           Merge reads with bbmerge before UMI extraction (default: false)
--umi_len          UMI length in bp (default: 12) - Omit this parameter to disable UMI extraction
--umi_regex        UMI regex pattern, expect UMI downstream of this seq (default: ATCGTCGGA)

Help:

```bash
nextflow run ./nextflow_pipeline/seqy.nf --help
```

## Troubleshooting

If you encounter any issues:

1. **Environment Setup Problems**: Ensure you have conda installed and the config script has execute permissions (`chmod +x config.sh`)
2. **Missing Dependencies**: Re-run the setup script or manually install missing packages with `conda install` or `pip install`

## Requirements

- Nextflow
- Conda package manager

For additional support or questions, please refer to the application documentation or contact the development team.