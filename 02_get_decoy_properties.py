from rdkit import Chem
from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.rdmolops import GetFormalCharge
from concurrent.futures import ProcessPoolExecutor
from rdkit.Chem.rdMolDescriptors import (CalcNumRotatableBonds,
                                         CalcNumHBA,
                                         CalcNumHBD)


def get_properties(smiles):
    # XXX: ExactMolWt includes H-atom weight. Should it?
    # XXX: CalcNumRotatableBonds is not "strict". Should it be?
    # XXX: DUD-E authors used a different program that calculated miLogP. Is our
    # proxy OK?

    try:
        mol = Chem.MolFromSmiles(smiles)
    except:
        return ""

    properties = [ExactMolWt(mol), MolLogP(mol),
                  CalcNumRotatableBonds(mol),
                  CalcNumHBA(mol), CalcNumHBD(mol),
                  GetFormalCharge(mol)]
    return ",".join([str(prop) for prop in properties])


def get_properties_batch(smiles_batch):
    return [f"{smiles},{get_properties(smiles)}" for smiles in smiles_batch]


print("Appending properties to data/all_zinc.smi")

with open("data/all_zinc.smi", "r") as f:
    all_zinc = [smiles.strip("\n") for smiles in f.readlines()]

indices = list(range(0, len(all_zinc), len(all_zinc) // 80)) + [-1]
input_subsets = []
for i in range(len(indices) - 1):
    if indices[i + 1] != -1:
        input_subsets.append(all_zinc[indices[i]: indices[i + 1]])
    else:
        input_subsets.append(all_zinc[indices[i]:])

with ProcessPoolExecutor() as executor:
    with_properties = executor.map(get_properties_batch, input_subsets)

with_properties = [line for sublist in with_properties
                   for line in sublist if line.split(",")[1]]

print(f"{len(all_zinc) - len(with_properties)} invalid compounds removed.")

with open("data/all_zinc_w_prop.csv", "w") as f:
    f.write("\n".join(with_properties))

print("Written data/all_zinc_w_prop.csv")

