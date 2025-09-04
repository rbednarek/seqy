#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * NGS Analysis Pipeline
 * Workflow: FastQC -> Optional Merge -> UMI Extract -> Alignment
 */

// Pipeline parameters
params.samplesheet = null
params.outdir = './results'
params.reference_genome = null
params.merged = false
params.umi_len = null
params.umi_regex = 'ATCGTCGGA'
params.help = false

// Help message
if (params.help) {
    log.info """
    NGS Analysis Pipeline
    ====================
    
    Usage:
    nextflow run pipeline.nf --samplesheet samples.csv --reference_genome genome.fa
    
    Required parameters:
    --samplesheet       CSV file with columns: sample_name,R1,R2
    --reference_genome  Path to reference genome FASTA file
    
    Optional parameters:
    --outdir           Output directory (default: ./results)
    --merged           Merge reads with bbmerge before UMI extraction (default: false)
    --umi_len          UMI length in bp (default: 12) - Omit this parameter to disable UMI extraction
    --umi_regex        UMI regex pattern, expect UMI seq downstream  (default: ATCGTCGGA) 
    --help             Show this help message
    """
    exit 0
}

// Validate required parameters
if (!params.samplesheet) {
    error "Please specify a samplesheet with --samplesheet"
}
if (!params.reference_genome) {
    error "Please specify a reference genome with --reference_genome"
}

/*
 * WORKFLOW
 */
workflow {
    
    // Read samplesheet and create channel
    samples_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            def sample_name = row.sample_name
            def r1 = file(row.R1)
            def r2 = file(row.R2)
            
            // Validate files exist
            if (!r1.exists()) error "R1 file does not exist: ${r1}"
            if (!r2.exists()) error "R2 file does not exist: ${r2}"
            
            return [sample_name, r1, r2]
        }
    
    // Reference genome
    reference_ch = Channel.fromPath(params.reference_genome)
    
    // Pipeline steps
    FASTQC(samples_ch)
    
    // Merge-Aware UMI extraction (optional - contingent on umi len param)
    if (params.merged) {
        BBMERGE(samples_ch)
        if (params.umi_len) {
            UMI_EXTRACT_MERGED(BBMERGE.out)
            processed_reads_ch = UMI_EXTRACT_MERGED.out.reads
        } else {
            processed_reads_ch = BBMERGE.out
        }
    } else {
        if (params.umi_len) {
            UMI_EXTRACT_PAIRED(samples_ch)
            processed_reads_ch = UMI_EXTRACT_PAIRED.out.reads
        } else {
            processed_reads_ch = samples_ch
        }
    }
    
    // Alignment (Bowtie2)
    ALIGN(processed_reads_ch, reference_ch)

    // UMI deduplication (optional - contingent on umi len param)
    if (params.umi_len) {
        if (params.merged) {
            UMI_DEDUP_MERGED(ALIGN.out)
        } else {
            UMI_DEDUP_PAIRED(ALIGN.out)
        }
    }
    
    // MultiQC
    MULTIQC(file("${params.outdir}"))
}

/*
 * PROCESSES
 */

process FASTQC {
    tag "$sample_name"
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    
    input:
    tuple val(sample_name), path(r1), path(r2)
    
    output:
    path "*.{html,zip}"
    
    script:
    """
    fastqc -t $task.cpus $r1 $r2
    """
}

process BBMERGE {
    tag "$sample_name"
    publishDir "${params.outdir}/merged", mode: 'copy'
    
    input:
    tuple val(sample_name), path(r1), path(r2)
    
    output:
    tuple val(sample_name), path("${sample_name}_merged.fastq.gz")
    
    script:
    """
    bbmerge.sh \\
        in1=$r1 \\
        in2=$r2 \\
        out=${sample_name}_merged.fastq.gz \\
        outu1=${sample_name}_unmerged_R1.fastq.gz \\
        outu2=${sample_name}_unmerged_R2.fastq.gz \\
        threads=$task.cpus
    """
}

process UMI_EXTRACT_MERGED {
    tag "$sample_name"
    publishDir "${params.outdir}/umi_extracted", mode: 'copy'
    
    input:
    tuple val(sample_name), path(merged_reads)
    
    output:
    tuple val(sample_name), path("${sample_name}_extracted.fastq.gz"), emit: reads
    path "${sample_name}_umi_extract.log", emit: log
    
    script:
    """
    umi_tools extract \\
        --extract-method=regex \\
        --bc-pattern=".{${params.umi_len}}(?=${params.umi_regex})" \\
        --stdin=$merged_reads \\
        --stdout=${sample_name}_extracted.fastq.gz \\
        2> ${sample_name}_umi_extract.log
    """
}

process UMI_EXTRACT_PAIRED {
    tag "$sample_name"
    publishDir "${params.outdir}/umi_extracted", mode: 'copy'
    
    input:
    tuple val(sample_name), path(r1), path(r2)
    
    output:
    tuple val(sample_name), path("${sample_name}_extracted_R1.fastq.gz"), path("${sample_name}_extracted_R2.fastq.gz"), emit: reads
    path "${sample_name}_umi_extract.log", emit: log
    
    script:
    """
    umi_tools extract \\
        --extract-method=regex \\
        --bc-pattern=".{${params.umi_len}}(?=${params.umi_regex})" \\
        --stdin=$r1 \\
        --stdout=${sample_name}_extracted_R1.fastq.gz \\
        --read2-in=$r2 \\
        --read2-out=${sample_name}_extracted_R2.fastq.gz \\
        2> ${sample_name}_umi_extract.log
    """
}

process ALIGN {
    tag "$sample_name"
    publishDir "${params.outdir}/alignments", mode: 'copy'
    
    input:
    tuple val(sample_name), path(reads)
    path reference
    
    output:
    tuple val(sample_name), path("${sample_name}.bam"), path("${sample_name}.bam.bai")
    
    script:
    if (params.merged) {
        """
        # Build Bowtie2 index if it doesn't exist
        if [ ! -f ${reference}.1.bt2 ]; then
            bowtie2-build $reference ${reference}
        fi
        
        # Align merged reads using Bowtie2
        bowtie2 -p $task.cpus \\
            -x $reference \\
            -U $reads | \\
            samtools view -bS - | \\
            samtools sort -@ $task.cpus -o ${sample_name}.bam -
        
        # Index BAM file
        samtools index ${sample_name}.bam
        
        # Generate alignment statistics
        samtools flagstat ${sample_name}.bam > ${sample_name}_flagstat.txt
        samtools stats ${sample_name}.bam > ${sample_name}_stats.txt
        """
    } else {
        """
        # Build Bowtie2 index if it doesn't exist
        if [ ! -f ${reference}.1.bt2 ]; then
            bowtie2-build $reference ${reference}
        fi
        
        # Align paired reads using Bowtie2
        bowtie2 -p $task.cpus \\
            -x $reference \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} | \\
            samtools view -bS - | \\
            samtools sort -@ $task.cpus -o ${sample_name}.bam -
        
        # Index BAM file
        samtools index ${sample_name}.bam
        
        # Generate alignment statistics
        samtools flagstat ${sample_name}.bam > ${sample_name}_flagstat.txt
        samtools stats ${sample_name}.bam > ${sample_name}_stats.txt
        """
    } 
}

process UMI_DEDUP_MERGED {
    tag "$sample_name"
    publishDir "${params.outdir}/dedup_bam", mode: 'copy'
    
    input:
    tuple val(sample_name), path("${sample_name}.bam"), path("${sample_name}.bam.bai")
    
    output:
    tuple val(sample_name), path("${sample_name}_dedup.bam"), path("${sample_name}_dedup.bam.bai")
    path "${sample_name}_dedup.log", emit: log
    
    script:
    """
    umi_tools dedup \\
        -I ${sample_name}.bam \\
        -S ${sample_name}_dedup.bam \\
        2> ${sample_name}_dedup.log

    samtools index ${sample_name}_dedup.bam    
    """
}

process UMI_DEDUP_PAIRED {
    tag "$sample_name"
    publishDir "${params.outdir}/dedup_bam", mode: 'copy'
    
    input:
    tuple val(sample_name), path("${sample_name}.bam"), path("${sample_name}.bam.bai")
    
    output:
    tuple val(sample_name), path("${sample_name}_dedup.bam"), path("${sample_name}_dedup.bam.bai")
    path "${sample_name}_dedup.log", emit: log
    
    script:
    """
    umi_tools dedup \\
        --paired \\
        -I ${sample_name}.bam \\
        -S ${sample_name}_dedup.bam \\
        2> ${sample_name}_dedup.log

    samtools index ${sample_name}_dedup.bam    
    """
}


process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    
    input:
    path fastqc_files
    
    output:
    path "multiqc_report.html"
    path "multiqc_data"
    
    script:
    """
    multiqc . --outdir . --force
    """
}

/*
 * WORKFLOW COMPLETION
 */
workflow.onComplete {
    log.info """
    Pipeline completed!
    Merged reads: ${params.merged}
    UMI extraction: ${params.umi_len}
    Results are in: ${params.outdir}
    """
}