# Snakefile
import pandas as pd

SAMPLES = ["patient_A", "patient_B", "patient_C", "patient_D", "patient_E", "patient_F"]

rule all:
    input:
        "results/cox_ph_results.csv",
        "results/km_plot.png"

rule index_ref:
    input: "data/ref.fasta"
    output: "data/ref.fasta.bwt"
    shell: "bwa index {input}"

rule bwa_mem:
    input:
        ref="data/ref.fasta",
        idx="data/ref.fasta.bwt", # Ensures indexing finishes first
        fq="data/{sample}.fastq"
    output: "mapped/{sample}.bam"
    shell: "bwa mem {input.ref} {input.fq} | samtools view -Sb - > {output}"

rule sort_bam:
    input: "mapped/{sample}.bam"
    output: "mapped/{sample}.sorted.bam"
    shell: "samtools sort {input} -o {output}"

rule index_bam:
    input: "mapped/{sample}.sorted.bam"
    output: "mapped/{sample}.sorted.bam.bai"
    shell: "samtools index {input}"

rule variant_calling:
    input:
        ref="data/ref.fasta",
        bam="mapped/{sample}.sorted.bam",
        bai="mapped/{sample}.sorted.bam.bai"
    output: "vcf/{sample}.vcf"
    shell: "bcftools mpileup -Ou -f {input.ref} {input.bam} | bcftools call -mv -Ov -o {output}"

rule extract_features:
    input: expand("vcf/{sample}.vcf", sample=SAMPLES)
    output: "results/genomic_matrix.csv"
    script: "scripts/extract_features.py"

rule survival_analysis:
    input:
        matrix="results/genomic_matrix.csv",
        clinical="data/clinical.csv"
    output:
        res="results/cox_ph_results.csv",
        plot="results/km_plot.png"
    script: "scripts/survival_analysis.py"