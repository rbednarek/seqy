# Seq-Y
`/'siː.ki/`

Pipeline written in Nextflow and Snakemake for standard NGS analysis. Includes read processing, alignment and optional UMI-based deduplication functionality.  

## Installation

### Step 1: Clone the Repository

```bash
git clone https://github.com/rbednarek/seqy.git
cd seqy
```


### Step 2: Set Up Conda Environment

Run the config script to create a new conda environment and install all required dependencies.

**Nextflow:**
```bash
bash ./config.sh nextflow
```
**Snakemake:**
```bash
bash ./config.sh snakemake
```
This script will:
- Create a new conda environment named `seqy`
- Install workflow manager and other required software packages
- **Note:** If you want to install both nextflow and snakemake in this environment simply run `bash ./config.sh`

### Step 3: Activate the Envrionment

After the setup script completes, activate the conda environment if necessary:

```bash
conda activate seqy
```

## Nextflow Pipeline

Basic Usage:

```bash
nextflow run ./nextflow/seqy.nf --samplesheet samples.csv --reference_genome genome.fa
```
    
Required parameters:
- **--samplesheet** - CSV file with columns: sample_name,R1,R2
- **--reference_genome** - Path to reference genome FASTA file

    
Optional parameters:
- **--outdir** - Output directory (default: ./results)
- **--merged** - Merge reads with bbmerge before UMI extraction (default: false)
- **--umi_len** - UMI length in bp (default: 12) - Omit this parameter to disable UMI extraction
- **--umi_regex** - UMI regex pattern, expect UMI downstream of this seq (default: ATCGTCGGA)

Help:

```bash
nextflow run ./nextflow/seqy.nf --help
```

## Snakemake Pipeline

This version requires a YAML input file in addition to a Samplesheet.

Basic Usage:

```bash
snakemake --snakefile ./snakemake/seqy.snakefile --configfile config.yaml --cores 4
```

Command Line Parameters:
- **--configfile** - Config file in `yaml` format
- **--cores** - Designate # of cores to use

Config File Parameters:
- **reference:** - Path to reference genome (ex. reference/genome.fa) **(Required)**
- **samplesheet:** - CSV file with columns: sample_name,R1,R2 **(Required)**
- **outdir:** - Results directory **(Required)**
- **merged:** - Merge reads with bbmerge before UMI extraction **(Optional)**
- **umi_len:** - UMI length in bp - Omit this parameter to disable UMI extraction **(Optional)**
- **umi_regex:** - UMI regex pattern, expect UMI downstream of this seq **(Optional)**

Help:
```bash
snakemake --help
```

An example config.yaml file can be found here: `./test_data/snk_config.yaml`).

## Example Dataset for Pipeline Testing
For testing purposes, several paired-end FASTQ files with UMIs have been generated alongside a reference genome. Files can be found here: `./test_data`

In the main `seqy` directory run one of the following:

**Nextflow (with UMI Deduplication):**:
```bash
nextflow run ./nextflow/seqy.nf --samplesheet ./test_data/samplesheet.csv --reference_genome ./test_data/genome/genome.fa --merged true --umi_len 12 --umi_regex ATCGTCGGA
```

**Snakemake (with UMI Deduplication):**
```bash
snakemake --snakefile  ./snakemake/seqy.snakefile --configfile ./test_data/snk_config.yaml --cores 4
```

## Troubleshooting

If you encounter any issues:

1. **Environment Setup Problems**: Ensure you have conda installed and the config script has execute permissions (`chmod +x config.sh`)
2. **Missing Dependencies**: Re-run the setup script or manually install missing packages with `conda install` or `pip install`

## Requirements

- Nextflow/Snakemake
- Conda package manager

For additional support or questions, please refer to the application documentation or contact the development team.


