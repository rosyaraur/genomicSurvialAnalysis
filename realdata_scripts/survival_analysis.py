import pandas as pd
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter, KaplanMeierFitter

# Load data
df_expr = pd.read_csv(snakemake.input.matrix)
df_clin = pd.read_csv(snakemake.input.clinical)
df = pd.merge(df_clin, df_expr, on="sample")

# Define the target gene you are investigating
# (For a real dataset, this would be an actual Ensembl or HGNC ID like 'ENSG00000141510' or 'TP53')
TARGET_GENE = "Gene_A" 

# 1. Fit Cox Proportional Hazards Model (using continuous TPM data)
cph = CoxPHFitter()
model_df = df[['time', 'status', TARGET_GENE, 'age']]
cph.fit(model_df, duration_col='time', event_col='status')
cph.summary.to_csv(snakemake.output.res)

# 2. Kaplan-Meier Plot (Stratified by Median Expression)
kmf = KaplanMeierFitter()
plt.figure(figsize=(8, 6))

# Calculate median to split the cohort
median_expr = df[TARGET_GENE].median()
high_mask = df[TARGET_GENE] > median_expr

# Plot High Expression
if high_mask.any():
    kmf.fit(df['time'][high_mask], event_observed=df['status'][high_mask], label="High Expression")
    kmf.plot_survival_function()

# Plot Low Expression
if (~high_mask).any():
    kmf.fit(df['time'][~high_mask], event_observed=df['status'][~high_mask], label="Low Expression")
    kmf.plot_survival_function()

plt.title(f"Survival Analysis Stratified by {TARGET_GENE} Expression")
plt.ylabel("Survival Probability")
plt.xlabel("Time (Months)")
plt.savefig(snakemake.output.plot)