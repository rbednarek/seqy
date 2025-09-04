import pandas as pd

config_yaml: "config.yaml"

samples = pd.read_csv(config["samplesheet"]).set_index("sample_name").T.to_dict()

if config.get('umi_len'):
    final_bam = expand("{outdir}/dedup/{sample}_dedup.bam", outdir=config["outdir"], sample=samples)
    final_bai = expand("{outdir}/dedup/{sample}_dedup.bam.bai", outdir=config["outdir"], sample=samples)
else:
    final_bam = expand("{outdir}/bam/{sample}.bam", outdir=config["outdir"], sample=samples)
    final_bai = expand("{outdir}/bam/{sample}.bam.bai", outdir=config["outdir"], sample=samples)

rule all:
    input:
        expand("{outdir}/multiqc/multiqc_report.html", outdir=config["outdir"]),
        final_bam,
        final_bai

if config.get('merged'):
    exec('use rule bbmerge as merge')


if config.get('umi_len'):
    if config.get('merged'):
        exec('use rule umi_extract_merged as umi_extract')
        exec('use rule align_merged as align')
        exec('use rule dedup as dedup')
    else:
        exec('use rule umi_extract_paired as umi_extract')
        exec('use rule align_paired as align')


rule fastqc:
    input:
        r1=lambda wildcards: samples[wildcards.sample]["R1"],
        r2=lambda wildcards: samples[wildcards.sample]["R2"]
    output:
        html1="{outdir}/fastqc/{sample}_R1_fastqc.html",
        html2="{outdir}/fastqc/{sample}_R2_fastqc.html"
    params:
        outdir="{outdir}/fastqc"
    threads: 2
    shell:
        """
        fastqc -t {threads} -o {params.outdir} {input.r1} {input.r2}
        """

rule bbmerge:
    input:
        r1=lambda wildcards: samples[wildcards.sample]["R1"],
        r2=lambda wildcards: samples[wildcards.sample]["R2"]
    output:
        merged="{outdir}/merged/{sample}_merged.fastq.gz",
        unmerged_r1="{outdir}/merged/{sample}_unmerged_R1.fastq.gz",
        unmerged_r2="{outdir}/merged/{sample}_unmerged_R2.fastq.gz"
    threads: 4
    shell:
        """
        bbmerge.sh in1={input.r1} in2={input.r2} \
            out={output.merged} \
            outu1={output.unmerged_r1} \
            outu2={output.unmerged_r2} \
            threads={threads}
        """

rule umi_extract_merged:
    input:
        merged="{outdir}/merged/{sample}_merged.fastq.gz"
    output:
        fastq="{outdir}/umi/{sample}_extracted.fastq.gz"
    params:
        pattern=config["umi_regex"],
        len=config["umi_len"]
    shell:
        """
        umi_tools extract \
            --extract-method=regex \
            --bc-pattern="(?<={params.pattern}).{{{params.len}}}" \
            --stdin={input.merged} \
            --stdout={output.fastq} \
            2> {outdir}/umi/{sample}_extracted.log
        """

rule umi_extract_paired:
    input:
        r1=lambda wildcards: samples[wildcards.sample]["R1"],
        r2=lambda wildcards: samples[wildcards.sample]["R2"]
    output:
        r1="{outdir}/umi/{sample}_R1_extracted.fastq.gz",
        r2="{outdir}/umi/{sample}_R2_extracted.fastq.gz"
    params:
        pattern=config["umi_regex"],
        len=config["umi_len"]
    shell:
        """
        umi_tools extract \
            --extract-method=regex \
            --bc-pattern="(?<={params.pattern}).{{{params.len}}}" \
            --stdin={input.r1} \
            --stdout={output.r1} \
            --read2-in={input.r2} \
            --read2-out={output.r2}
        """

rule align_merged:
    input:
    fq = lambda wc: (f"{config['outdir']}/umi/{wc.sample}_extracted.fastq.gz"
        if config.get("umi_len")
        else f"{config['outdir']}/merged/{wc.sample}_merged.fastq.gz"),        
    ref=config["reference"]
    output:
        bam="{outdir}/bam/{sample}.bam",
        bai="{outdir}/bam/{sample}.bam.bai"
    threads: 4
    shell:
        """
        if [ ! -f {input.ref}.1.bt2 ]; then
            bowtie2-build {input.ref} {input.ref}
        fi
        bowtie2 -p {threads} -x {input.ref} -U {input.fq} | samtools view -bS - |
            samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

rule align_paired:
    input:
        r1=lambda wc: (
            f"{config['outdir']}/umi/{wc.sample}_R1_extracted.fastq.gz"
            if config.get("umi_len")
            else samples[wc.sample]["R1"]
        ),
        r2=lambda wc: (
            f"{config['outdir']}/umi/{wc.sample}_R2_extracted.fastq.gz"
            if config.get("umi_len")
            else samples[wc.sample]["R2"]
        ),
        ref=config["reference"]
    output:
        bam="{outdir}/bam/{sample}.bam",
        bai="{outdir}/bam/{sample}.bam.bai"
    threads: 4
    shell:
        """
        if [ ! -f {input.ref}.1.bt2 ]; then
            bowtie2-build {input.ref} {input.ref}
        fi
        bowtie2 -p {threads} -x {input.ref} -1 {input.r1} -2 {input.r2} | samtools view -bS - |
            samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

rule dedup:
    input:
        bam="{outdir}/bam/{sample}.bam",
        bai="{outdir}/bam/{sample}.bam.bai"
    output:
        bam="{outdir}/dedup/{sample}_dedup.bam",
        bai="{outdir}/dedup/{sample}_dedup.bam.bai"
    params:
        paired = "--paired" if not config.get('merged') else ""
    shell:
        """
        umi_tools dedup \
            {params.paired} \
            -I {input.bam} \
            -S {output.bam} \
            2> {outdir}/dedup/{sample}_dedup.log
        samtools index {output.bam}
        """

rule multiqc:
    input:
        expand("{outdir}/fastqc/{sample}_R1_fastqc.html", outdir=config["outdir"], sample=samples)
    output:
        html="{outdir}/multiqc/multiqc_report.html"
    shell:
        """
        multiqc {config['outdir']} -o {config['outdir']}/multiqc
        """
