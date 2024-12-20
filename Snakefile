from pathlib import Path

samples_dir = Path(config["samples_dir"])
genome_dir = Path(config["genome_dir"])
results_dir = Path(config["results_dir"])


def genome(wildcards):
    return [
        genome_dir / x
        for x in [
            "genome.fa",
            "genome.fa.amb",
            "genome.fa.ann",
            "genome.fa.bwt",
            "genome.fa.fai",
            "genome.fa.pac",
            "genome.fa.sa",
        ]
    ]


SAMPLES = ["A", "B"]


rule all:
    input:
        results_dir / "plots/quals.svg"

rule bwa_map:
    input:
        genome=genome,
        fastq=samples_dir / "{sample}.fastq"
    output:
        results_dir / "mapped_reads/{sample}.bam"
    shell:
        "bwa mem {input.genome[0]} {input.fastq} | samtools view -Sb - > {output}"


rule samtools_sort:
    input:
        results_dir / "mapped_reads/{sample}.bam"
    output:
        results_dir / "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} -O bam {input} > {output}"


rule samtools_index:
    input:
        results_dir /  "sorted_reads/{sample}.bam"
    output:
        results_dir /  "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"


rule bcftools_call:
    input:
        genome=genome,
        bam=expand(results_dir / "sorted_reads/{sample}.bam", sample=SAMPLES),
        bai=expand(results_dir / "sorted_reads/{sample}.bam.bai", sample=SAMPLES)
    output:
        results_dir / "calls/all.vcf"
    shell:
        "bcftools mpileup -f {input.genome[0]} {input.bam} | "
        "bcftools call -mv - > {output}"


rule plot_quals:
    input:
        results_dir / "calls/all.vcf"
    output:
        results_dir / "plots/quals.svg"
    script:
        "scripts/plot-quals.py"
