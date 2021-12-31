"""Lorem ipsum."""

# Steps:
# 1. Generate 50 decoys for each ligand in the *expanded* ligands dictionary
# through property-matching to get between 3000 and 9000 candidate decoys.
# 2. Across the candidates for each ligand, compute the maximum Tc (based on
# ECFP4 fingerprints) between the candidate and any one associated ligand.
# 3. Keep only the top 25% most dissimilar candidates.
# 4. Deduplication thing.

#need to have shuffling

import os
import argparse
from shutil import rmtree
import pickle
import multiprocessing

import numpy as np
import pandas as pd
import SharedArray as sa

from rdkit import Chem, DataStructs
from progressbar import progressbar
from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem.rdmolops import GetFormalCharge
from concurrent.futures import ProcessPoolExecutor
from rdkit.Chem.rdMolDescriptors import (CalcExactMolWt,
                                         CalcNumRotatableBonds,
                                         CalcNumHBA,
                                         CalcNumHBD)



np.random.seed(12345)

rmtree("data/sa", ignore_errors=True)
os.makedirs("data/sa")
if not os.path.isdir("data/cache"):
    os.makedirs("data/cache")

COLUMNS = ("smiles",
           "mol_wt", "logp", "num_rotatable", "num_hba", "num_hbd", "net_charge")
MWT_RANGES =  [ 20,  35,  50,  65,  80, 100, 125]
LOGP_RANGES = [0.4, 0.8, 1.2, 1.8, 2.4, 3.0, 3.6]
RB_RANGES =   [  1,   2,   2,   3,   3,   4,   5]
NHA_RANGES =  [  0,   1,   2,   2,   3,   3,   4]
NHD_RANGES =  [  0,   0,   1,   1,   2,   2,   3]
CHG_RANGES =  [  0,   0,   0,   0,   0,   1,   2]
RANGES = [MWT_RANGES, LOGP_RANGES, RB_RANGES, NHA_RANGES, NHD_RANGES, CHG_RANGES]

BATCH_HELP_MESSAGE = \
"""The size of the subset of input ligands to be processed concurrently. If the
total number of ligands submitted for decoy candidate generation is 160, for
example, a batch size of 80 would create two batches. The ligands in the first
batch would have candidates generated for them concurrently, with the results
cached, before proceeding to the second batch.

Depending on the number of CPU cores in your system, it is recommended not to
make the batch size too small to make this program's parallelization
counterproductive due to its overhead: the default size is thus 128. However, a
smaller batch size may be preferred when this program has to be stopped/started
repeatedly, since batches that have been completed are saved and skipped on
subsequent runs."""


def create_sa(prop):
    smiles = np.load(f"data/sorted/{prop}_smiles_reorder.npy")
    values = np.load(f"data/sorted/{prop}_values.npy")

    smiles_ = sa.create(f"file://data/sa/{prop}_smiles", smiles.shape, dtype=smiles.dtype)
    values_ = sa.create(f"file://data/sa/{prop}_values", values.shape, dtype=values.dtype)

    np.copyto(smiles_, smiles)
    np.copyto(values_, values)

    return smiles_, values_


def get_batched_ligand_pairs(input_path, batch_size):
    """Split referenced ligand pairs list into batches of batch_size."""

    with open(input_path, "r") as f:
       ligand_pairs = [line.strip("\n") for line in f.readlines()]

    batched_ligand_pairs = []
    indices = list(range(0, len(ligand_pairs), batch_size)) + [-1]
    for i in range(len(indices) - 1):
        if indices[i + 1] == -1:
            batched_ligand_pairs.append(ligand_pairs[indices[i]:indices[i + 1]])
        else:
            batched_ligand_pairs.append(ligand_pairs[indices[i]:])

    return batched_ligand_pairs



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


def get_properties(mol):
    """Get the six to-match properties for the Mol object."""

    return np.array([CalcExactMolWt(mol), MolLogP(mol),
                     CalcNumRotatableBonds(mol, strict=True),
                     CalcNumHBA(mol), CalcNumHBD(mol),
                     GetFormalCharge(mol)])


def get_candidates(smiles):
    """Get candidate decoy indices for a single ligand."""

    mol = Chem.MolFromSmiles(smiles)
    ligand_props = get_properties(mol)

    prop_matched_indices = []
    window_index = 0
    # XXX: What role does min have here? I think we can put a break if we've hit
    # at least 3000.
    while len(prop_matched_indices) < 9000 and window_index < 7:
        latest = get_property_matched(ligand_props,
                                      window_index)
        window_index += 1

        if latest is None:
            continue

        if not len(prop_matched_indices):
            tried_union = latest.copy()
        else:
            tried_union = np.union1d(prop_matched_indices, latest)

        if len(tried_union) > 9000:
            latest = np.setdiff1d(latest, prop_matched_indices,
                                  assume_unique=True)
            try:
                prop_matched_indices = \
                    np.union1d(prop_matched_indices,
                               np.random.choice(latest,
                                                9000 - len(prop_matched_indices),
                                                replace=False))
            except Exception as e:
                #print(e)
                prop_matched_indices = np.union1d(latest,
                                                   prop_matched_indices)
        else:
            prop_matched_indices = tried_union.copy()

        if len(tried_union) >= 3000:
            break

    print(len(prop_matched_indices))

    return prop_matched_indices


def write_candidate_dict(batched_ligand_pairs, num_cores, job_name):
    """Write the candidate dictionary for the dataset."""
    batches_to_do = [batch_number
                     for batch_number in range(len(batched_ligand_pairs))
                     if not os.path.isfile(f"data/cache/{job_name}/{batch_number}")]

    prefix = "Resuming" if batches_to_do[0] else "Starting"
    print(f"{prefix} job {job_name} from batch {batches_to_do[0]}/{len(batches_to_do)}...")

    # TODO: remove serial alt when parallel confirmed.
    #"""
    candidate_dict = dict()
    #"""
    for curr_batch in progressbar(batches_to_do):
        #"""
        for smiles_pair in progressbar(batched_ligand_pairs[curr_batch]):
            candidate_dict[smiles_pair] = \
                get_candidates(smiles_pair.split("_")[1])
        #"""
        """
        with ProcessPoolExecutor(max_workers=num_cores) as executor:
            candidates = executor.map(get_candidates,
                                      batched_ligand_pairs)

        curr_candidate_dict = dict()
        curr_smiles = list(LIGAND_DICT[pdb_id].keys())
        for i, subset in enumerate(candidates):
            curr_candidate_dict[curr_smiles[i]] = subset
        """

        with open(f"data/candidates/{pdb_id}_candidate_dict.pkl", "wb") as f:
            pickle.dump(curr_candidate_dict, f)

    print("Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Supply input and save paths.")
    # TODO: remove default after testing and make required=True.
    parser.add_argument("job_name", default="test",
                        help=(f"Job name to be used to set the results cache "
                              f"location."))
    # TODO: remove default after testing and make required=True.
    parser.add_argument("-i",
                        default="../get_data/BindingDB/bindingdb_actives_protonated",
                        help=(f"Relative path to input actives_protonated "
                              f"SMILES file as created by 04_get_protonated.py"))
    parser.add_argument("--batch_size", default=128, help=BATCH_HELP_MESSAGE)
    parser.add_argument("--num_cores", default=multiprocessing.cpu_count(),
                        help=(f"Number of CPU cores to use for processing every"
                              f" batch Default: all available cores."))
    args = vars(parser.parse_args())

    if not os.path.isdir(f"data/cache/{args['job_name']}"):
        os.makedirs(f"data/cache/{args['job_name']}")

    print(f"Creating decoy candidates for {args['i']}")

    SORTED_ZINC_DFS = [create_sa(prop) for prop in COLUMNS[1:]]

    print("[1/2] Loading data...")
    batched_ligand_pairs = get_batched_ligand_pairs(args["i"], args["batch_size"])

    print(f"[2/2] Getting decoy candidates...")
    write_candidate_dict(batched_ligand_pairs,
                         args["num_cores"],
                         args["job_name"])

    print("Done.")

