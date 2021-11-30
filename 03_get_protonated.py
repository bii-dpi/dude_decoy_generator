import pickle
import dimorphite_dl

import numpy as np

from rdkit import Chem, RDLogger
from collections import defaultdict
from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem.Descriptors import ExactMolWt
from concurrent.futures import ProcessPoolExecutor
from rdkit.Chem.rdMolDescriptors import (CalcNumRotatableBonds,
                                         CalcNumHBA,
                                         CalcNumHBD)


RDLogger.DisableLog("rdApp.*")

PATHS = ["../get_data/BindingDB/bindingdb_actives",
         "../get_data/DUDE/dude_actives"]


def get_mol_dict(path):
    print("[1/3] Getting Mol object for each ligand...")
    with open(path, "r") as f:
        smiles = [line.split()[:2] for line in f.readlines()]
    mols = [[Chem.MolFromSmiles(pair[0]), pair[0], pair[1]]
            for pair in smiles]

    mol_dict = defaultdict(list)
    for i in range(len(mols)):
        mol_dict[mols[i][-1]].append(mols[i][:2])

    return mol_dict


def get_properties(mol):
    # XXX: ExactMolWt includes H-atom weight. Should it?
    # XXX: CalcNumRotatableBonds is not "strict". Should it be?
    # XXX: DUD-E authors used a different program that calculated miLogP. Is our
    # proxy OK?
    # TODO: Their sixth descriptor was "added net charge". What does that mean?

    return np.array([ExactMolWt(mol), MolLogP(mol),
                     CalcNumRotatableBonds(mol),
                     CalcNumHBA(mol), CalcNumHBD(mol)])


def get_protonated_single(mol):
    protonated_mols = dimorphite_dl.run_with_mol_list(
        [mol],
        min_ph=6.0,
        max_ph=8.0,
        silent=True
    )

    orig_properties = get_properties(mol)
    included_mols = [mol]
    included_properties = [orig_properties]
    for other_mol in protonated_mols:
        curr_properties = get_properties(other_mol)
        add = True
        for included_properties_ in included_properties:
            if (included_properties_ == curr_properties).all():
                add = False
                break

        if add:
            included_mols.append(other_mol)

    return [Chem.MolToSmiles(other_mol)
            for other_mol in included_mols]


def get_protonated_dict_protein(curr_mols):
    ligand_dict = dict()
    for mol, smiles in curr_mols:
        ligand_dict[smiles] = get_protonated_single(mol)
    return ligand_dict


def save_protonated_dict(path):
    print(f"Saving dict for {path}")
    mol_dict = get_mol_dict(path)
    protonated_dict = dict()
    print("[2/3] Protonating ligands for each protein...")
    with ProcessPoolExecutor() as executor:
        protonated_subdicts = executor.map(get_protonated_dict_protein,
                                           [mol_dict[pdb_id]
                                            for pdb_id in mol_dict.keys()])
    protonated_dict = dict(zip(mol_dict.keys(), protonated_subdicts))

    print("[3/3] Saving...")
    save_path = f"data/{path.split('/')[-1]}_protonated_dict.pkl"
    with open(save_path, "wb") as f:
        pickle.dump(protonated_dict, f)

    print(f"Saved in {save_path}.\n")


for path in PATHS:
    save_protonated_dict(path)

