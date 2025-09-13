setup-all: nextflow snakemake

nextflow:
	bash setup/config.sh nextflow

snakemake:
	bash setup/config.sh snakemake	

pre-commit:
	pre-commit install


	
