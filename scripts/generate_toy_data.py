# generate_toy_data.py
import os
import random
import pandas as pd

os.makedirs("data", exist_ok=True)

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

    # Toy FASTQ: introduce random variations so bcftools finds mutations
    with open(f"data/{s}.fastq", "w") as f:
        for i in range(50):
            start = random.randint(0, 900)
            read = list(ref_seq[start:start+50])
            if random.random() > 0.9: 
                read[25] = random.choice("ACGT") # Mutate a base
            f.write(f"@{s}_{i}\n{''.join(read)}\n+\n{'I'*50}\n")

pd.DataFrame(clin_data).to_csv("data/clinical.csv", index=False)
print("Toy data generated in ./data/")
