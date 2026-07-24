# scripts/extract_features.py
import pandas as pd
import os

vcf_files = snakemake.input
output_file = snakemake.output[0]

data = []
for vcf in vcf_files:
    sample_name = os.path.basename(vcf).replace(".vcf", "")
    
    # Parse VCF to count variants (ignoring header lines)
    with open(vcf) as f:
        var_count = sum(1 for line in f if not line.startswith("#"))
    
    # Binary feature: 1 if patient has mutations, 0 if wildtype
    has_mut = 1 if var_count > 0 else 0
    data.append({"sample": sample_name, "has_mutation": has_mut, "total_variants": var_count})

df = pd.DataFrame(data)
df.to_csv(output_file, index=False)