# Genomic survival analysis pipeline 

1. **Configure the Conda Environment:** Installs all required bioinformatics and data science tools.
Create a file named `environment.yaml` in your root directory. This includes all the tools we discovered were necessary during troubleshooting.

```yaml
name: survival_ngs
channels:
  - conda-forge
  - bioconda
dependencies:
  - snakemake-minimal
  - bwa
  - samtools
  - bcftools
  - python=3.10
  - pandas
  - lifelines
  - matplotlib

```

Build and activate the environment:

```bash
conda env create -f environment.yaml
conda activate survival_ngs

```


2. **Create the Directory Structure:**
Set up your project folders to ensure Snakemake knows where to find inputs and place outputs. Run this in your terminal:

```bash
mkdir -p data mapped vcf results scripts

```


3. **Simulate the Toy Data:** Includes the adjusted mutation rate for statistical variance.
Create `scripts/generate_toy_data.py`. This script generates a tiny reference genome, clinical metadata for 6 patients, and toy FASTQ files. The mutation threshold is set to `0.9` to ensure we get a mix of mutated and wildtype patients.

```python
import os
import random
import pandas as pd

# 1. Generate Toy Reference Genome
ref_seq = "".join(random.choices("ACGT", k=1000))
with open("data/ref.fasta", "w") as f:
    f.write(">chr1\n" + ref_seq + "\n")

# 2. Generate Clinical Data and FASTQ files
samples = ["patient_A", "patient_B", "patient_C", "patient_D", "patient_E", "patient_F"]
clin_data = []

for s in samples:
    # Random survival data (time in months, status 1=event/0=censored)
    time = random.randint(10, 100)
    status = random.choice([0, 1])
    age = random.randint(40, 80)
    clin_data.append({"sample": s, "time": time, "status": status, "age": age})

    # Toy FASTQ: 10% chance to mutate a base (0.9 threshold)
    with open(f"data/{s}.fastq", "w") as f:
        for i in range(50):
            start = random.randint(0, 900)
            read = list(ref_seq[start:start+50])
            if random.random() > 0.9: 
                read[25] = random.choice("ACGT")
            f.write(f"@{s}_{i}\n{''.join(read)}\n+\n{'I'*50}\n")

pd.DataFrame(clin_data).to_csv("data/clinical.csv", index=False)
print("Toy data generated in ./data/")

```


4. **Define the Snakemake Pipeline:**
Create a file named `Snakefile` (no extension) in your root directory. This orchestrates the flow from raw reads to the final survival plot.

```python
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
        idx="data/ref.fasta.bwt",
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

```


5. **Write the Feature Extraction Script:**
Create `scripts/extract_features.py`. This script converts the VCF files into a binary matrix (1 if mutations exist, 0 if wildtype).

```python
import pandas as pd
import os

vcf_files = snakemake.input
output_file = snakemake.output[0]

data = []
for vcf in vcf_files:
    sample_name = os.path.basename(vcf).replace(".vcf", "")
    
    with open(vcf) as f:
        var_count = sum(1 for line in f if not line.startswith("#"))
    
    has_mut = 1 if var_count > 0 else 0
    data.append({"sample": sample_name, "has_mutation": has_mut, "total_variants": var_count})

df = pd.DataFrame(data)
df.to_csv(output_file, index=False)

```


6. **Write the Survival Analysis Script:** Includes the zero-variance safeguard.
Create `scripts/survival_analysis.py`. This script merges the clinical and genomic data, checks for statistical viability, and generates the Cox regression and Kaplan-Meier plots.

```python
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter, KaplanMeierFitter
import sys

# Load data
df_gen = pd.read_csv(snakemake.input.matrix)
df_clin = pd.read_csv(snakemake.input.clinical)
df = pd.merge(df_clin, df_gen, on="sample")

# Variance check
if df['has_mutation'].nunique() <= 1:
    print("WARNING: 'has_mutation' has zero variance (all samples are the same).", file=sys.stderr)
    
    pd.DataFrame(columns=["covariate", "exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]).to_csv(snakemake.output.res, index=False)
    
    plt.figure(figsize=(8, 6))
    plt.text(0.5, 0.5, 'Not enough variance for Kaplan-Meier\n(All samples have same mutation status)', 
             horizontalalignment='center', verticalalignment='center')
    plt.savefig(snakemake.output.plot)
    sys.exit(0)

# 1. Fit Cox Proportional Hazards Model
cph = CoxPHFitter()
model_df = df[['time', 'status', 'has_mutation', 'age']]
cph.fit(model_df, duration_col='time', event_col='status')
cph.summary.to_csv(snakemake.output.res)

# 2. Plot Kaplan-Meier Curves
kmf = KaplanMeierFitter()
plt.figure(figsize=(8, 6))

mask = df['has_mutation'] == 1
if mask.any():
    kmf.fit(df['time'][mask], event_observed=df['status'][mask], label="Mutated")
    kmf.plot_survival_function()

if (~mask).any():
    kmf.fit(df['time'][~mask], event_observed=df['status'][~mask], label="Wildtype")
    kmf.plot_survival_function()

plt.title("Survival Analysis by Mutation Status")
plt.ylabel("Survival Probability")
plt.xlabel("Time (Months)")
plt.savefig(snakemake.output.plot)

```


7. **Execute the Pipeline:**
Run the entire workflow from your project's root directory:

```bash
# 1. Generate the raw files
python scripts/generate_toy_data.py

# 2. Run the pipeline
snakemake --cores 1

```