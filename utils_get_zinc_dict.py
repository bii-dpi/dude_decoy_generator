"""Combine .smi files in data/raw_files/ in a random order."""

import os
from pickle import dump
from progressbar import progressbar


print("Combining all downloaded .smi files into dictionary...")


def read_smi(fname):
    """Get SMILES from fname."""

    assert fname.endswith(".smi"), f"{fname} not valid .smi file."

    with open(f"data/raw_files/{fname}", "r") as f:
        return [line.strip("\n").split() for line in f.readlines()[1:]]


# Compile list of all unique ZINC15 SMILES.
all_smiles = []
for fname in progressbar(os.listdir("data/raw_files")):
    all_smiles += read_smi(fname)
all_smiles = dict(all_smiles)

with open("data/zinc_dict.pkl", "wb") as f:
    dump(all_smiles, f)

print(f"Written ZINC SMILES-ID dictionary.")

