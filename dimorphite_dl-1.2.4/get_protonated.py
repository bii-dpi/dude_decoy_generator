import dimorphite_dl

from tqdm import tqdm
from rdkit import Chem


PATHS = ["../get_data/BindingDB/bindingdb_actives",
         "../get_data/DUDE/dude_actives"]


def get_mols(path):
    with open(path, "r") as f:
        smiles = [line.split()[:2] for line in f.readlines()]
    return [[Chem.MolFromSmiles(pair[0]), pair[1]]
            for pair in smiles]

protonated_mols = dimorphite_dl.run_with_mol_list(
    get_mols(PATHS[1])[:10],
    min_ph=6.0,
    max_ph=8.0,
)
print([Chem.MolToSmiles(m) for m in protonated_mols])

