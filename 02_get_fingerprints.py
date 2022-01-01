"""Write the ECFP4 fingerprint for each SMILES in all_zinc.smi"""

import os
from pickle import dump
from rdkit.Chem import AllChem
from rdkit import Chem, DataStructs
from concurrent.futures import ProcessPoolExecutor


print("Getting ECFP4 fingerprints for all ZINC15 compounds...")


if not os.path.isdir("data/fingerprints"):
    os.makedirs("data/fingerprints")
if not os.path.isdir("data/fingerprints/ecfp4"):
    os.makedirs("data/fingerprints/ecfp4")


def get_fingerprint(smiles):
    """Get the ECFP4 fingerprint for the input SMILES."""

    # If the input SMILES is invalid, return None.
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception as e:
        print(f"{smiles} could not be made to Mol: {e}")
        return

    return AllChem.GetMorganFingerprint(mol, 2)


def get_fingerprints_batch(pair):
    """Get ECFP4 fingerprints for each SMILES in the input batch."""

    smiles_batch, i = pair

    with open(f"data/fingerprints/ecfp4/{i}.pkl", "wb") as f:
        dump([get_fingerprint(smiles) for smiles in smiles_batch], f)


with open("data/all_zinc.smi", "r") as f:
    all_zinc = [smiles.strip("\n") for smiles in f.readlines()]

# Break the complete SMILES list into batches to be processed in parallel.
indices = list(range(0, len(all_zinc), 10000)) + [-1]
input_subsets = []
for i in range(len(indices) - 1):
    if indices[i + 1] != -1:
        input_subsets.append((all_zinc[indices[i]:indices[i + 1]], indices[i]))
    else:
        input_subsets.append((all_zinc[indices[i]:], indices[i]))


if __name__ == "__main__":
    with ProcessPoolExecutor() as executor:
        executor.map(get_fingerprints_batch, input_subsets)

print("Done.")

