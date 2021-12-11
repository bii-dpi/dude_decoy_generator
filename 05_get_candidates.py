# Steps:
# 1. Generate 50 decoys for each ligand in the *expanded* ligands dictionary
# through property-matching to get between 3000 and 9000 candidate decoys.
# 2. Across the candidates for each ligand, compute the maximum Tc (based on
# ECFP4 fingerprints) between the candidate and any one associated ligand.
# 3. Keep only the top 25% most dissimilar candidates.
# 4. Deduplication thing.

#need to have shuffling

import os
import shutil
import pickle
import multiprocessing

import numpy as np
import pandas as pd
import SharedArray as sa

from rdkit import Chem, DataStructs
from progressbar import progressbar
from collections import defaultdict
from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem.Descriptors import ExactMolWt
from concurrent.futures import ProcessPoolExecutor
from rdkit.Chem.rdMolDescriptors import (CalcNumRotatableBonds,
                                         CalcNumHBA,
                                         CalcNumHBD)



np.random.seed(12345)

shutil.rmtree("data/sa", ignore_errors=True)
os.makedirs("data/sa")


DATASET = "DUDE"


def create_sa(prop):
    smiles = np.load(f"data/sorted/{prop}_smiles_reorder.npy")
    values = np.load(f"data/sorted/{prop}_values.npy")

    smiles_ = sa.create(f"file://data/sa/{prop}_smiles", smiles.shape, dtype=smiles.dtype)
    values_ = sa.create(f"file://data/sa/{prop}_values", values.shape, dtype=values.dtype)

    np.copyto(smiles_, smiles)
    np.copyto(values_, values)

    return smiles_, values_


def get_smiles_in_range(value, delta, curr_df):
    left_pos = np.searchsorted(curr_df[1], value - delta, side="left")
    right_pos = np.searchsorted(curr_df[1], value + delta, side="right")

    if left_pos == len(curr_df[1]):
        left_pos -= 1
    if right_pos == len(curr_df[1]):
        right_pos -= 1

    assert right_pos >= left_pos

    return curr_df[0][left_pos: right_pos + 1]


def get_property_matched(curr_ligand_props, window_index):
    all_matched = None
    for i in range(len(COLUMNS[1:])):
        if all_matched is None:
            all_matched = get_smiles_in_range(curr_ligand_props[i],
                                              RANGES[i][window_index],
                                              SORTED_ZINC_DFS[i])
        else:
            all_matched = np.intersect1d(all_matched,
                                         get_smiles_in_range(curr_ligand_props[i],
                                                             RANGES[i][window_index],
                                                             SORTED_ZINC_DFS[i]),
                                         assume_unique=True)
            if len(all_matched) == 0:
                return None
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
    if len(smiles_list) > 10:
        smiles_list = [smiles_list[0]] + \
                        np.random.choice(smiles_list[1:],
                                         size=9, replace=False).tolist()

    ligand_mols = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]
    ligand_props = [get_mol_properties(mol) for mol in ligand_mols]

    all_property_matched = []
    for curr_ligand_props in ligand_props:
        window_index = 0
        curr_property_matched = []
        while len(curr_property_matched) < 50 and window_index < 7:
            latest = get_property_matched(curr_ligand_props,
                                          window_index)
            window_index += 1

            if latest is None:
                continue

            if not len(curr_property_matched):
                tried_union = latest.copy()
            else:
                tried_union = np.union1d(curr_property_matched, latest)

            if len(tried_union) > 50:
                latest = np.setdiff1d(latest, curr_property_matched,
                                      assume_unique=True)
                try:
                    curr_property_matched = \
                        np.union1d(curr_property_matched,
                                   np.random.choice(latest,
                                                    50 - len(curr_property_matched),
                                                    replace=False))
                except Exception as e:
                    #print(e)
                    curr_property_matched = np.union1d(latest,
                                                       curr_property_matched)
            else:
                curr_property_matched = tried_union.copy()
        if not len(all_property_matched):
            all_property_matched = curr_property_matched.copy()
        else:
            all_property_matched = np.union1d(curr_property_matched,
                                              all_property_matched)
        #print("prop", len(all_property_matched))

    #print("prop", len(all_property_matched))

    return all_property_matched


def write_candidate_dict(all_pdb_ids):
    """Write the candidate dictionary for the dataset."""
    all_pdb_ids = [pdb_id for pdb_id in all_pdb_ids
                   if not os.path.isfile(f"data/candidates/{pdb_id}_candidate_dict.pkl")]
    """
    candidate_dict = defaultdict(dict)
    """
    for pdb_id in progressbar(all_pdb_ids):
        """
        for smiles in progressbar(LIGAND_DICT[pdb_id]):
            candidate_dict[pdb_id][smiles] = \
                get_candidates(LIGAND_DICT[pdb_id][smiles])
        """
        with ProcessPoolExecutor(max_workers=80) as executor:
            candidates = executor.map(get_candidates,
                                      [LIGAND_DICT[pdb_id][smiles]
                                       for smiles in LIGAND_DICT[pdb_id]])

        curr_candidate_dict = dict()
        curr_smiles = list(LIGAND_DICT[pdb_id].keys())
        for i, subset in enumerate(candidates):
            curr_candidate_dict[curr_smiles[i]] = subset
        with open(f"data/candidates/{pdb_id}_candidate_dict.pkl", "wb") as f:
            pickle.dump(curr_candidate_dict, f)

    print("Done.")


if __name__ == "__main__":
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
    print("[1/2] Loading data...")
    LIGAND_DICT = pd.read_pickle(f"data/{DATASET.lower()}_actives_protonated_dict.pkl")
    SORTED_ZINC_DFS = [create_sa(prop) for prop in COLUMNS[1:]]

    print(f"[2/2] Getting candidates for {len(LIGAND_DICT)} proteins...")
    write_candidate_dict(LIGAND_DICT.keys())

