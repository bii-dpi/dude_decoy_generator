# Steps:
# 1. Generate 50 decoys for each ligand in the *expanded* ligands dictionary
# through property-matching to get between 3000 and 9000 candidate decoys.
# 2. Across the candidates for each ligand, compute the maximum Tc (based on
# ECFP4 fingerprints) between the candidate and any one associated ligand.
# 3. Keep only the top 25% most dissimilar candidates.
# 4. Deduplication thing.

#need to have shuffling


import numpy as np
import pandas as pd

from bisect import bisect
from rdkit import Chem, DataStructs


DATASET = "BindingDB"

ATD_RATIO = 50
COLUMNS = ("smiles",
           "mol_wt", "logp", "num_rotatable", "num_hba", "num_hbd", "net_charge")
MWT_RANGES =  [ 20,  35,  50,  65,  80, 100, 125]
LOGP_RANGES = [0.4, 0.8, 1.2, 1.8, 2.4, 3.0, 3.6]
RB_RANGES =   [  1,   2,   2,   3,   3,   4,   5]
NHA_RANGES =  [  0,   1,   2,   2,   3,   3,   4]
NHD_RANGES =  [  0,   0,   1,   1,   2,   2,   3]
CHG_RANGES =  [  0,   0,   0,   0,   0,   1,   2]
RANGES = [MWT_RANGES, LOGP_RANGES, RB_RANGES, NHA_RANGES, NHD_RANGES, CHG_RANGES]


print(f"Creating decoys for {DATASET}.")
print("[1/X] Loading data...")
LIGAND_DICT = pd.read_pickle(f"data/{DATASET.lower()}_actives_protonated_dict.pkl")
SORTED_ZINC_DFS = [pd.read_csv(f"data/all_zinc_w_prop_{prop}.csv")
                   for prop in COLUMNS[1:]]
PROP_MEDIANS = [df.iloc[1].median() for df in SORTED_ZINC_DFS]

MAX_CAND_TC = dict(zip(SORTED_ZINC_DFS[0].smiles,
                       [-1 for _ in SORTED_ZINC_DFS[0].smiles]))


def get_properties(mol):
    return np.array([ExactMolWt(mol), MolLogP(mol),
              CalcNumRotatableBonds(mol),
              CalcNumHBA(mol), CalcNumHBD(mol),
              get_net_charge(Chem.MolToSmiles(mol))])


def get_indices_in_range(value, delta, sorted_arr, median):
    selected_indices = []
    min_value = value - delta
    max_value = value + delta
    if value < median:
        for i, other_value in enumerate(sorted_arr):
            if other_value > max_value:
                break
            elif other_value < min_value:
                continue
            selected_indices.append(i)
    else:
        for i in range(len(sorted_arr) - 1, -1, -1):
            other_value = sorted_arr[i]
            if other_value > max_value:
                continue
            elif other_value < min_value:
                break
            selected_indices.append(i)

    return set(selected_indices)


def get_indices_in_range(value, delta, curr_df, median):
    sorted_arr = curr_df.iloc[:, 1].tolist()
    left_pos = bisect_left(sorted_arr, value - delta)
    right_pos = bisect_right(sorted_arr, value + delta)

    while sorted_arr[left_pos] < value - delta:
        left_pos += 1
    while sorted_arr[right_pos] > value + delta:
        right_pos -= 1

    return set(curr_df.iloc[range(left_pos, right_pos + 1)])


def get_property_matched_single(i, prop_value, window_index):
    get_indices_in_range(prop_value,
                         RANGES[i][window_index],
                         SORTED_ZINC_DFS[i][COLUMNS[i + 1]],
                         PROP_MEDIANS[i]),


def get_property_matched(ligand_mol, ligand_props, window_index):
    all_matched = set(SORTED_ZINC_DFS[0].smiles)
    for i in range(len(COLUMNS[1:])):
        all_matched &= get_property_matched_single(i, ligand_props[i], window_index)
    return all_matched


def get_decoys(pdb_id, ligand_smiles):
    ligand_mols = list(map(Chem.MolFromSmiles,
                           LIGAND_DICT[pdb_id][ligand_smiles]))
    ligand_props = [get_properties(mol) for mol in ligand_mols]
    """
    ligand_fingerprints = [AllChem.GetMorganFingerprint(ligand_mol, 2)
                           for ligand_mol in ligand_mols]
    """

    all_property_matched = set()
    window_index = 0
    while len(all_property_matched) <= 9000:
        for i, ligand_mol in enumerate(ligand_mols):
            all_property_matched |= get_property_matched(ligand_mol,
                                                         ligand_props[i],
                                                         window_index)
        window_index += 1
    return all_property_matched


def get_sim(active_fingerprint, candidate_smiles):
    """Get Tc between the compounds' ECP4 fingerprints."""
    candidate_fingerprint = AllChem.GetMorganFingerprint(
                                Chem.MolFromSmiles(candidate_smiles),
                                2)
    return DataStructs.DiceSimilarity(active_fingerprint,
                                      candidate_fingerprint)


def get_curr_candidates(all_candidates, active_fingerprint):
    """Get the bottom 25% of candidates by Tc to the active."""
    all_sims = [get_sim(active_fingerprint, candidate_smiles)
                for candidate_smiles in all_candidates]
    return all_candidates[np.argsort(all_sims)][:int(len(all_candidates) * 0.25)]


def get_decoys(all_candidates, active_smiles):
    """Get 50 decoys for single SMILES."""
    print("Getting dissimilar-enough candidates...")
    active_fingerprint = AllChem.GetMorganFingerprint(
                            Chem.MolFromSmiles(active_smiles),
                            2)
    curr_candidates = get_curr_candidates(all_candidates, active_fingerprint)

    assert len(curr_candidates) >= 50, \
            f"{active_smiles} doesn't have enough dissimilar candidates."




