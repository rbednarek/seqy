import pandas as pd

configfile: "config.yaml"

samples = pd.read_csv(config["samplesheet"]).set_index("sample_name").T.to_dict()

if config.get('umi_len'):
    final_bam = expand("{outdir}/dedup/{sample}_dedup.bam", outdir=config["outdir"], sample=samples)
    final_bai = expand("{outdir}/dedup/{sample}_dedup.bam.bai", outdir=config["outdir"], sample=samples)
else:
    if config.get('merged'):
        final_bam = expand("{outdir}/bam_merged/{sample}.bam", outdir=config["outdir"], sample=samples)
        final_bai = expand("{outdir}/bam_merged/{sample}.bam.bai", outdir=config["outdir"], sample=samples)
    else:
        final_bam = expand("{outdir}/bam_paired/{sample}.bam", outdir=config["outdir"], sample=samples)
        final_bai = expand("{outdir}/bam_paired/{sample}.bam.bai", outdir=config["outdir"], sample=samples)

rule all:
    input:
        expand("{outdir}/multiqc/multiqc_report.html", outdir=config["outdir"]),
        final_bam,
        final_bai

# Rule order to resolve ambiguity
ruleorder: umi_extract_merged > umi_extract_paired
ruleorder: align_merged > align_paired

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
        fastq="{outdir}/umi/{sample}_extracted.fastq.gz",
        log="{outdir}/umi/{sample}_extracted.log"
    params:
        pattern=config.get("umi_regex"),
        len=config.get("umi_len")
    shell:
        """
        umi_tools extract \
            --extract-method=regex \
            --bc-pattern="(?P<umi_1>.{{{params.len}}})(?P<discard_1>{params.pattern})" \
            --stdin={input.merged} \
            --stdout={output.fastq} \
            2> {output.log}
        """

rule umi_extract_paired:
    input:
        r1=lambda wildcards: samples[wildcards.sample]["R1"],
        r2=lambda wildcards: samples[wildcards.sample]["R2"]
    output:
        r1="{outdir}/umi/{sample}_R1_extracted.fastq.gz",
        r2="{outdir}/umi/{sample}_R2_extracted.fastq.gz",
        log="{outdir}/umi/{sample}_extracted.log"
    params:
        pattern=config.get("umi_regex"),
        len=config.get("umi_len")
    shell:
        """
        umi_tools extract \
            --extract-method=regex \
            --bc-pattern="(?P<umi_1>.{{{params.len}}})(?P<discard_1>{params.pattern})" \
            --stdin={input.r1} \
            --stdout={output.r1} \
            --read2-in={input.r2} \
            --read2-out={output.r2} \
            2> {output.log}
        """

rule align_merged:
    input:
        fq=lambda wc: (f"{config['outdir']}/umi/{wc.sample}_extracted.fastq.gz"
            if config.get("umi_len")
            else f"{config['outdir']}/merged/{wc.sample}_merged.fastq.gz"),        
        ref=config["reference"]
    output:
        bam="{outdir}/bam_merged/{sample}.bam",  # Different directory
        bai="{outdir}/bam_merged/{sample}.bam.bai"
    threads: 4
    shell:
        """
        if [ ! -f {input.ref}.1.bt2 ]; then
            bowtie2-build {input.ref} {input.ref}
        fi
        bowtie2 -p {threads} --local -x {input.ref} -U {input.fq} | samtools view -bS - | \
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
        bam="{outdir}/bam_paired/{sample}.bam",  # Different directory
        bai="{outdir}/bam_paired/{sample}.bam.bai"
    threads: 4
    shell:
        """
        if [ ! -f {input.ref}.1.bt2 ]; then
            bowtie2-build {input.ref} {input.ref}
        fi
        bowtie2 -p {threads} --local -x {input.ref} -1 {input.r1} -2 {input.r2} | samtools view -bS - | \
            samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

rule dedup:
    input:
        bam=lambda wc: (f"{config['outdir']}/bam_merged/{wc.sample}.bam"
                        if config.get('merged')
                        else f"{config['outdir']}/bam_paired/{wc.sample}.bam"),
        bai=lambda wc: (f"{config['outdir']}/bam_merged/{wc.sample}.bam.bai"
                        if config.get('merged')
                        else f"{config['outdir']}/bam_paired/{wc.sample}.bam.bai")
    output:
        bam="{outdir}/dedup/{sample}_dedup.bam",
        bai="{outdir}/dedup/{sample}_dedup.bam.bai",
        log="{outdir}/dedup/{sample}_dedup.log"
    params:
        paired="--paired" if not config.get('merged') else ""
    shell:
        """
        umi_tools dedup \
            {params.paired} \
            -I {input.bam} \
            -S {output.bam} \
            2> {output.log}
        samtools index {output.bam}
        """

rule multiqc:
    input:
        expand("{outdir}/fastqc/{sample}_R1_fastqc.html", outdir=config["outdir"], sample=samples)
    output:
        html="{outdir}/multiqc/multiqc_report.html"
    shell:
        """
        multiqc {config[outdir]} -o {config[outdir]}/multiqc
        """