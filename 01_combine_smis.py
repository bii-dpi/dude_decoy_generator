"""Combine .smi files in data/raw_files/ in a random order."""

import os
import numpy as np
from progressbar import progressbar


print("Combining all downloaded .smi files...")

np.random.seed(12345)


def read_smi(fname):
    """Get SMILES from fname."""

    assert fname.endswith(".smi"), f"{fname} not valid .smi file."

    with open(f"data/raw_files/{fname}", "r") as f:
        return [line.split()[0] for line in f.readlines()[1:]]


# Compile list of all unique ZINC15 SMILES.
all_smiles = []
for fname in progressbar(os.listdir("data/raw_files")):
    all_smiles += read_smi(fname)
all_smiles = np.unique(all_smiles)

np.random.shuffle(all_smiles)

with open("data/all_zinc.smi", "w") as f:
    f.write("\n".join(all_smiles))

print(f"Written {len(all_smiles)} ZINC SMILES.")

