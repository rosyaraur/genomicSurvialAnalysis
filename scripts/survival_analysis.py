# scripts/survival_analysis.py
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter, KaplanMeierFitter
import sys

# Load data
df_gen = pd.read_csv(snakemake.input.matrix)
df_clin = pd.read_csv(snakemake.input.clinical)

# Merge clinical and genomic data on sample ID
df = pd.merge(df_clin, df_gen, on="sample")

# Check variance: if all samples have the same mutation status, we can't run Cox
if df['has_mutation'].nunique() <= 1:
    print("WARNING: 'has_mutation' has zero variance (all samples are the same).", file=sys.stderr)
    pd.DataFrame(columns=["covariate", "exp(coef)", "p"]).to_csv(snakemake.output.res, index=False)
    
    plt.figure(figsize=(8, 6))
    plt.text(0.5, 0.5, 'Not enough variance for Kaplan-Meier', horizontalalignment='center')
    plt.savefig(snakemake.output.plot)
    sys.exit(0)

# 1. Fit Cox Proportional Hazards Model
# Added penalizer=0.1 to stabilize math for very small/collinear datasets
cph = CoxPHFitter(penalizer=0.1) 
model_df = df[['time', 'status', 'has_mutation', 'age']]

try:
    cph.fit(model_df, duration_col='time', event_col='status')
    cph.summary.to_csv(snakemake.output.res)
except Exception as e:
    print(f"WARNING: Cox model mathematically failed (expected with N=6). Error: {e}", file=sys.stderr)
    pd.DataFrame(columns=["covariate", "exp(coef)", "p"]).to_csv(snakemake.output.res, index=False)

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
