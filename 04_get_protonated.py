"""Write protonated forms of each active in supplied actives file."""

import argparse
import dimorphite_dl
import numpy as np
from rdkit import Chem
from collections import defaultdict
from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem.rdmolops import GetFormalCharge
from concurrent.futures import ProcessPoolExecutor
from rdkit.Chem.rdMolDescriptors import (CalcExactMolWt,
                                         CalcNumRotatableBonds,
                                         CalcNumHBA,
                                         CalcNumHBD)


np.random.seed(12345)


print("Saving the SMILES of each input active's protonated states for pH 6-8...")


def get_mols(input_path):
    """Get Mol objects of all actives SMILES in referenced file."""

    def get_mol(smiles):
        try:
            return Chem.MolFromSmiles(smiles)
        except Exception as e:
            print(f"{smiles} could not be processed: {e}")
            return

    with open(input_path, "r") as f:
        smiles = [line.strip("\n") for line in f.readlines()]

    mols = [[get_mol(curr_smiles), curr_smiles] for curr_smiles in smiles]
    old_len = len(mols)
    # Filter out SMILES that could not be converted to Mol objects.
    mols = [pair for pair in mols if pair[0]]
    if len(mols) != old_len:
        print(f"Warning: {old_len - len(mols)} actives SMILES could not be"
              f" processed.")

    return mols


def group_into_batches(mols):
    """Group the Mol objects in the input list into batches of 128."""

    batched_mols = []
    indices = list(range(0, len(mols), 128)) + [-1]
    for i in range(len(indices) - 1):
        if indices[i + 1] == -1:
            batched_mols.append(mols[indices[i]:])
        else:
            batched_mols.append(mols[indices[i]:indices[i + 1]])

    return batched_mols


def get_properties(mol):
    """Get the six properties associated with the input Mol."""

    return np.array([CalcExactMolWt(mol), MolLogP(mol),
                     CalcNumRotatableBonds(mol, strict=True),
                     CalcNumHBA(mol), CalcNumHBD(mol),
                     GetFormalCharge(mol)])


def get_protonated_single(mol, orig_smiles):
    """Get SMILES of the given Mol"s protonated states with unique property sets."""

    protonated_mols = dimorphite_dl.run_with_mol_list(
        [mol],
        min_ph=6.0,
        max_ph=8.0,
        silent=True
    )
    np.random.shuffle(protonated_mols)

    # Maintain a list of the Mol objects to keep for decoy generation and their
    # (unique) property-sets.
    included_mols = [mol]
    included_properties = [get_properties(mol)]
    for other_mol in protonated_mols:
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
    return [f"{orig_smiles}_{Chem.MolToSmiles(mol_)}"
            for mol_ in included_mols]


def get_protonated(curr_mol_batch):
    """Get SMILES pairs for protonated forms of actives in the input list."""

    smiles_pairs = []
    for mol, smiles in curr_mol_batch:
        smiles_pairs += get_protonated_single(mol, smiles)

    return smiles_pairs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Supply input and save paths.")
    # TODO: remove default after testing and make required=True.
    parser.add_argument("-i", default="data/test_actives",
                        help="Relative path to input actives SMILES file.")
    parser.add_argument("-o", default=None,
                        help=(f"Relative path to output actives and protonated "
                              f"states SMILES file. Default: input path + _protonated"))
    args = vars(parser.parse_args())

    input_path = args["i"]
    save_path = args["o"]
    if not save_path:
        save_path = input_path + "_protonated"
        print(f"No save path supplied; save path set to default: {save_path}")

    print("[1/3] Getting Mol object for each active...")
    mols = get_mols(input_path)

    print("[2/3] Protonating actives and saving those with unique property sets...")
    # Batch the mols into batches of 128 to justify parallelization overhead.
    batched_mols = group_into_batches(mols)

    with ProcessPoolExecutor() as executor:
        all_pairs = executor.map(get_protonated, batched_mols)
    # Flatten active-protonated SMILES pair list.
    all_pairs = [pair for sublist in all_pairs for pair in sublist]

    print(f"[3/3] Saving {len(all_pairs)} SMILES pairs created from "
          f"{len(mols)} actives...")
    with open(save_path, "w") as f:
        f.write("\n".join(all_pairs))

    print("Done.")

