"""Write relevant physicochemical properties for each SMILES in all_zinc.smi"""

from rdkit import Chem
from multiprocessing import cpu_count
from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem.rdmolops import GetFormalCharge
from concurrent.futures import ProcessPoolExecutor
from rdkit.Chem.rdMolDescriptors import (CalcExactMolWt,
                                         CalcNumRotatableBonds,
                                         CalcNumHBA,
                                         CalcNumHBD)


NUM_CORES = cpu_count()

print("Getting relevant physicochemical properties for all ZINC15 compounds...")


def get_properties(smiles):
    """Get a comma-separated line of the input SMILES compound's properties."""

    # If the input SMILES is invalid, return an empty string.
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception as e:
        print(f"{smiles} could not be made to Mol: {e}")
        return ""

    properties = [CalcExactMolWt(mol), MolLogP(mol),
                  CalcNumRotatableBonds(mol, strict=True),
                  CalcNumHBA(mol), CalcNumHBD(mol),
                  GetFormalCharge(mol)]
    return ",".join([str(prop) for prop in properties])


def get_properties_batch(smiles_batch):
    """Get property-lines of each SMILES in the input batch."""

    return [f"{smiles},{get_properties(smiles)}" for smiles in smiles_batch]


print("[1/1] Appending properties to data/all_zinc.smi ...")

with open("data/all_zinc.smi", "r") as f:
    all_zinc = [smiles.strip("\n") for smiles in f.readlines()]

# Break the complete SMILES list into batches to be processed in parallel.
indices = list(range(0, len(all_zinc), len(all_zinc) // NUM_CORES)) + [-1]
input_subsets = []
for i in range(len(indices) - 1):
    if indices[i + 1] != -1:
        input_subsets.append(all_zinc[indices[i]: indices[i + 1]])
    else:
        input_subsets.append(all_zinc[indices[i]:])

with ProcessPoolExecutor() as executor:
    with_properties = executor.map(get_properties_batch, input_subsets)

# Flatten the batched sublists in with_properties for valid lines.
with_properties = [line for sublist in with_properties
                   for line in sublist if line.split(",")[1]]

print(f"{len(all_zinc) - len(with_properties)} invalid compounds removed.")

with open("data/all_zinc_w_prop.csv", "w") as f:
    f.write("\n".join(with_properties))

print("Written data/all_zinc_w_prop.csv")

