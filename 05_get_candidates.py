# Steps:
# 1. Generate 50 decoys for each ligand in the *expanded* ligands dictionary
# through property-matching to get between 3000 and 9000 candidate decoys.
# 2. Across the candidates for each ligand, compute the maximum Tc (based on
# ECFP4 fingerprints) between the candidate and any one associated ligand.
# 3. Keep only the top 25% most dissimilar candidates.
# 4. Deduplication thing.

#need to have shuffling


import pickle

import numpy as np
import pandas as pd

from rdkit import Chem, DataStructs
from progressbar import progressbar
from collections import defaultdict
from rdkit.Chem.Crippen import MolLogP
from bisect import bisect_left, bisect_right
from rdkit.Chem.Descriptors import ExactMolWt
from concurrent.futures import ProcessPoolExecutor
from rdkit.Chem.rdMolDescriptors import (CalcNumRotatableBonds,
                                         CalcNumHBA,
                                         CalcNumHBD)

np.random.seed(12345)

DATASET = "DUDE"


def read_csv(prop):
    curr_df = pd.read_csv(f"data/all_zinc_w_prop_{prop}.csv")
    if prop in COLUMNS[3:]:
        curr_df[[prop]] = curr_df[[prop]].astype(int)

    return curr_df


def get_smiles_in_range(value, delta, curr_df):
    sorted_arr = curr_df.iloc[:, 1].tolist()
    left_pos = bisect_left(sorted_arr, value - delta)
    right_pos = bisect_right(sorted_arr, value + delta)

    while sorted_arr[left_pos] < value - delta:
        left_pos += 1
    while sorted_arr[right_pos] > value + delta:
        right_pos -= 1

    return set(curr_df.iloc[range(left_pos, right_pos + 1), 0].tolist())


def get_property_matched(curr_ligand_props, window_index):
    all_matched = None
    for i in range(len(COLUMNS[1:])):
        if all_matched is None:
            all_matched = get_smiles_in_range(curr_ligand_props[i],
                                              RANGES[i][window_index],
                                              SORTED_ZINC_DFS[i])
        else:
            all_matched &= get_smiles_in_range(curr_ligand_props[i],
                                               RANGES[i][window_index],
                                               SORTED_ZINC_DFS[i])

    return all_matched


def get_mol_properties(mol):
    """Get the six to-match properties for the Mol object."""
    def get_net_charge(mol):
        """Get molecule's net charge."""
        smiles = Chem.MolToSmiles(mol)
        pos = smiles.count("[+") + smiles.count("+]")
        neg = smiles.count("[-") + smiles.count("-]")

        return pos - neg

    return [ExactMolWt(mol), MolLogP(mol),
            CalcNumRotatableBonds(mol),
            CalcNumHBA(mol), CalcNumHBD(mol),
            get_net_charge(mol)]


def get_candidates(smiles_list):
    """Get candidate decoys for a single ligand."""
    ligand_mols = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]
    ligand_props = [get_mol_properties(mol) for mol in ligand_mols]

    all_property_matched = set()
    window_index = 0
    for curr_ligand_props in ligand_props:
        window_index = 0
        curr_property_matched = set()
        while len(curr_property_matched) < 50 and window_index < 7:
            latest = get_property_matched(curr_ligand_props,
                                          window_index)
            tried_union = curr_property_matched | latest
            if len(tried_union) > 50:
                latest -= curr_property_matched
                try:
                    curr_property_matched |= \
                        set(np.random.choice(list(latest),
                                             50 - len(curr_property_matched),
                                             replace=False).tolist())
                except:
                    curr_property_matched |= latest
            else:
                curr_property_matched = tried_union
            window_index += 1
        all_property_matched |= curr_property_matched
        print(len(all_property_matched))

    return all_property_matched


def write_candidate_dict(all_pdb_ids):
    """Write the candidate dictionary for the dataset."""
    candidate_dict = defaultdict(dict)
    for pdb_id in progressbar(all_pdb_ids):
        for smiles in progressbar(LIGAND_DICT[pdb_id]):
            candidate_dict[pdb_id][smiles] = \
                get_candidates(LIGAND_DICT[pdb_id][smiles])

    print("[3/3] Writing candidate dictionary...")
    with open(f"data/{DATASET.lower()}_candidate_dict.pkl", "wb") as f:
        pickle.dump(candidate_dict, f)

    print("Done.")


COLUMNS = ("smiles",
           "mol_wt", "logp", "num_rotatable", "num_hba", "num_hbd", "net_charge")
MWT_RANGES =  [ 20,  35,  50,  65,  80, 100, 125]
LOGP_RANGES = [0.4, 0.8, 1.2, 1.8, 2.4, 3.0, 3.6]
RB_RANGES =   [  1,   2,   2,   3,   3,   4,   5]
NHA_RANGES =  [  0,   1,   2,   2,   3,   3,   4]
NHD_RANGES =  [  0,   0,   1,   1,   2,   2,   3]
CHG_RANGES =  [  0,   0,   0,   0,   0,   1,   2]
RANGES = [MWT_RANGES, LOGP_RANGES, RB_RANGES, NHA_RANGES, NHD_RANGES, CHG_RANGES]


print(f"Creating decoy candidates for {DATASET}.")
print("[1/3] Loading data...")
LIGAND_DICT = pd.read_pickle(f"data/{DATASET.lower()}_actives_protonated_dict.pkl")
SORTED_ZINC_DFS = [read_csv(prop) for prop in COLUMNS[1:]]

print(f"[2/3] Getting candidates for {len(LIGAND_DICT)} proteins...")
write_candidate_dict(LIGAND_DICT.keys())

