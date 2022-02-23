"""
Write protonated forms of each active in input file which have a unique property
set.
"""

import argparse
import numpy as np
from rdkit import Chem
from collections import defaultdict
from progressbar import progressbar
from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem.rdmolops import GetFormalCharge
from rdkit.Chem.rdMolDescriptors import (CalcExactMolWt,
                                         CalcNumRotatableBonds,
                                         CalcNumHBA,
                                         CalcNumHBD)


np.random.seed(12345)


print("Saving qualifying SMILES of each input active's protonated states for pH 6-8...")


def get_mol(smiles):
    """Get Mol object from input SMILES."""

    try:
        return Chem.MolFromSmiles(smiles)
    except Exception as e:
        print(f"{smiles} could not be processed: {e}")
        return


def get_actives_dict(input_path):
    """Get Mol objects of all actives SMILES in referenced file."""

    with open(input_path, "r") as f:
        data = [line.strip("\n").split() for line in f.readlines()]

    mols_dict = {b_id: [get_mol(smiles), smiles] for smiles, b_id in data}
    old_len = len(mols_dict)
    # Filter out SMILES that could not be converted to Mol objects.
    mols_dict = {b_id: pair for b_id, pair in mols_dict.items() if pair[0]}
    if len(mols_dict) != old_len:
        print(f"Warning: {old_len - len(mols)} actives SMILES could not be"
              f" processed.")

    return mols_dict


def get_protonated_actives_dict(input_path):
    """Get Mol objects of all ligand SMILES in referenced file."""

    with open(input_path, "r") as f:
        data = [line.strip("\n").split() for line in f.readlines()]

    mols_dict = defaultdict(list)
    for smiles, b_id in data:
        curr_mol = get_mol(smiles)
        if curr_mol is not None:
            mols_dict[b_id].append([curr_mol, smiles])

    for b_id in mols_dict:
        np.random.shuffle(mols_dict[b_id])

    return mols_dict


def get_properties(mol):
    """Get the six properties associated with the input Mol."""

    return np.array([CalcExactMolWt(mol), MolLogP(mol),
                     CalcNumRotatableBonds(mol, strict=True),
                     CalcNumHBA(mol), CalcNumHBD(mol),
                     GetFormalCharge(mol)])


def select_protonated(active_mol, active_smiles, protonated_list):
    """Get SMILES of the given Mol"s protonated states with unique property sets."""

    # Maintain a list of the Mol objects to keep for decoy generation and their
    # (unique) property-sets.
    included_mols = [active_mol]
    included_properties = [get_properties(active_mol)]
    for other_mol, _ in protonated_list:
        curr_properties = get_properties(other_mol)
        add = True
        # Keep only Mol objects that have a sextuplet not seen in the
        # current properties list.
        for included_properties_ in included_properties:
            if (included_properties_ == curr_properties).all():
                add = False
                break
        if add:
            included_mols.append(other_mol)
            included_properties.append(curr_properties.copy())

    # Prepend the original active's SMILES to that of its protonated form. This
    # will later help in identifying which active decoy candidates correspond
    # to.
    return [f"{active_smiles}_{Chem.MolToSmiles(mol_)}"
            for mol_ in included_mols]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Supply input and save paths.")
    parser.add_argument("--input_actives", "-a", required=True,
                        help="Relative path to input actives SMILES file.")
    parser.add_argument("--input_actives_protonated", "-p", required=True,
                        help="Relative path to input protonated actives SMILES file.")
    parser.add_argument("--output_actives_protonated", "-o", default=None,
                        help=(f"Relative path to output actives and protonated "
                              f"states SMILES file. Default: input path + "
                              f"_protonated_selected"))
    args = vars(parser.parse_args())

    actives_path = args["input_actives"]
    protonated_actives_path = args["input_actives_protonated"]
    save_path = args["output_actives_protonated"]
    if not save_path:
        save_path = actives_path + "_protonated_selected"
        print(f"No save path supplied; save path set to default: {save_path}")

    print("[1/3] Getting Mol object for each active and its protonated forms...")
    actives_dict = get_actives_dict(actives_path)
    protonated_actives_dict = get_protonated_actives_dict(protonated_actives_path)

    print("[2/3] Selecting protonated actives with unique property sets...")
    all_pairs = []
    for b_id in progressbar(actives_dict):
        all_pairs += select_protonated(*actives_dict[b_id],
                                       protonated_actives_dict[b_id])

    print(f"[3/3] Saving {len(all_pairs)} SMILES pairs created from "
          f"{len(actives_dict)} actives...")
    with open(save_path, "w") as f:
        f.write("\n".join(all_pairs))

    print("Done.")

