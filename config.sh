#!/bin/bash

# Miniconda and seqy environment setup script
# This script installs miniconda, creates the seqy environment with required packages,
# and configures automatic activation

set -e

workflow=${1:-both}

echo "Starting miniconda and seqy environment setup..."


detect_os() {
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        echo "linux"
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        echo "macos"
    else
        echo "unsupported"
        exit 1
    fi
}


OS=$(detect_os)
echo "Detected OS: $OS"


if [[ "$OS" == "linux" ]]; then
    MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
elif [[ "$OS" == "macos" ]]; then
    if [[ $(uname -m) == "arm64" ]]; then
        MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh"
    else
        MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh"
    fi
fi


MINICONDA_DIR="$HOME/miniconda3"


if [[ -d "$MINICONDA_DIR" ]]; then
    echo "Miniconda already installed at $MINICONDA_DIR"
else
    echo "Installing miniconda..."

    wget -O /tmp/miniconda.sh "$MINICONDA_URL"
    
    bash /tmp/miniconda.sh -b -p "$MINICONDA_DIR"
    
    rm /tmp/miniconda.sh
    
    echo "Miniconda installed successfully"
fi

echo "Initializing conda..."
source "$MINICONDA_DIR/bin/activate"
conda init bash

source ~/.bashrc 2>/dev/null || source ~/.bash_profile 2>/dev/null || true

echo "Updating conda..."
conda update -n base -c defaults conda -y


base_pkgs=(
    python=3.10
    bowtie2
    bbmap
    umi_tools
    samtools
    fastqc
    multiqc
)

case "$workflow" in
    snakemake)
        workflow_pkgs=(snakemake)
        ;;
    nextflow)
        workflow_pkgs=(nextflow)
        ;;
    both|"")
        workflow_pkgs=(snakemake nextflow)
        ;;
    *)
        echo "Invalid option: $workflow"
        echo "Usage: $0 [snakemake|nextflow|both]"
        exit 1
        ;;
esac

conda create -n seqy -y -c conda-forge -c bioconda \
    "${base_pkgs[@]}" \
    "${workflow_pkgs[@]}"

echo "✅ Base environment created. Activating environment..."
conda activate seqy

echo "Activating seqy environment..."
source "$MINICONDA_DIR/bin/activate" seqy

echo "Configuring automatic activation of seqy environment..."

echo ""
echo "============================================================"
echo "Setup completed successfully!"
echo "============================================================"
echo ""
echo "Summary of what was installed:"
echo "• Miniconda3 at: $MINICONDA_DIR"
echo "• Environment: seqy"
echo "• Workflow: ${workflow}"
echo "• System packages: bowtie2, bbmap, umi_tools, samtools, fastqc, multiqc"
echo ""
echo "To manually activate the environment: conda activate seqy"
echo "To deactivate: conda deactivate"
echo ""
echo "Please restart your terminal or run: source ~/.bashrc"
echo "============================================================"