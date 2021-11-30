# Combine .smi files in data/raw_files/ randomly.
import os

import numpy as np

from tqdm import tqdm


np.random.seed(12345)


def read_smi(fname):
    assert fname.endswith(".smi"), f"{fname} not valid .smi file."

    with open(f"data/raw_files/{fname}", "r") as f:
        return [line.split()[0] for line in f.readlines()[1:]]


all_smiles = []
for fname in tqdm(os.listdir("data/raw_files")):
    all_smiles += read_smi(fname)
all_smiles = np.unique(all_smiles)

np.random.shuffle(all_smiles)

with open("data/all_zinc.smi", "w") as f:
    f.write("\n".join(all_smiles))

print(f"Written {len(all_smiles)} ZINC SMILES.")

