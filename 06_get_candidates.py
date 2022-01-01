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
import multiprocessing
import numpy as np
import pandas as pd
import SharedArray as sa
from pickle import dump
from shutil import rmtree
from collections import defaultdict
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
# TODO: Need to make the saved indices ints.

# Manage directories used for this program's storage.
rmtree("data/sa", ignore_errors=True)
os.makedirs("data/sa")

if not os.path.isdir("data/cache"):
    os.makedirs("data/cache")

if not os.path.isdir("data/candidate_dicts"):
    os.makedirs("data/candidate_dicts")

# Constants relating to properties decoy candidates are matched on.
COLUMNS = ("smiles",
           "mol_wt", "logp", "num_rotatable", "num_hba", "num_hbd", "net_charge")
MWT_RANGES =  [ 20,  35,  50,  65,  80, 100, 125]
LOGP_RANGES = [0.4, 0.8, 1.2, 1.8, 2.4, 3.0, 3.6]
RB_RANGES =   [  1,   2,   2,   3,   3,   4,   5]
NHA_RANGES =  [  0,   1,   2,   2,   3,   3,   4]
NHD_RANGES =  [  0,   0,   1,   1,   2,   2,   3]
CHG_RANGES =  [  0,   0,   0,   0,   0,   1,   2]
RANGES = [MWT_RANGES, LOGP_RANGES, RB_RANGES, NHA_RANGES, NHD_RANGES, CHG_RANGES]

# Help messages associated with command line arguments for this program.
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

RESET_HELP_MESSAGE = \
"""Reset the cache by deleting the stored results corresponding to the supplied
job name. For correct results, this flag must be raised if new data is being
supplied to the program or a different batch size is being used (under the same
job name)."""


def create_sa(prop):
    """Lorem ipsum."""
    smiles = np.load(f"data/sorted/{prop}_smiles_reorder.npy")
    values = np.load(f"data/sorted/{prop}_values.npy")

    smiles_ = sa.create(f"file://data/sa/{prop}_smiles", smiles.shape, dtype=smiles.dtype)
    values_ = sa.create(f"file://data/sa/{prop}_values", values.shape, dtype=values.dtype)

    np.copyto(smiles_, smiles)
    np.copyto(values_, values)

    return smiles_, values_


def get_batched_list(list_, batch_size):
    """Get the batched input list."""

    batched_list = []
    indices = list(range(0, len(list_), batch_size)) + [-1]
    for i in range(len(indices) - 1):
        if indices[i + 1] == -1:
            batched_list.append(list_[indices[i]:])
        else:
            batched_list.append(list_[indices[i]:indices[i + 1]])

    return batched_list


def get_batched_ligand_pairs(input_path, batch_size):
    """Get referenced ligand pairs list split into batches."""

    with open(input_path, "r") as f:
       ligand_pairs = [line.strip("\n") for line in f.readlines()]

    return get_batched_list(ligand_pairs, batch_size)


def get_smiles_in_range(value, delta, curr_df):
    """Lorem ipsum."""

    left_pos = np.searchsorted(curr_df[1], value - delta, side="left")
    right_pos = np.searchsorted(curr_df[1], value + delta, side="right")

    if left_pos == len(curr_df[1]):
        left_pos -= 1
    if right_pos == len(curr_df[1]):
        right_pos -= 1

    assert right_pos >= left_pos

    return curr_df[0][left_pos: right_pos + 1]


def get_property_matched(curr_ligand_props, window_index):
    """Lorem ipsum."""

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


def get_candidates_single(smiles_pair):
    """Get candidate decoy indices for a single ligand."""

    smiles = smiles_pair.split("_")[1]
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

#    print(len(prop_matched_indices))

    return prop_matched_indices


def get_candidates(smiles_pair_list):
    """Lorem ipsum."""

    """
    return [get_candidates_single(smiles_pair)
            for smiles_pair in smiles_pair_list]
    """
    return get_candidates_single(smiles_pair_list)


def write_candidate_dict(batched_ligand_pairs, num_cores, job_name):
    """Write the candidate dictionary for the dataset."""
    try:
        batched_candidate_dict = \
            pd.read_pickle(f"data/cache/{job_name}/batched_candidate_dict.pkl")
        batches_to_do = [i for i in range(len(batched_ligand_pairs))
                         if i not in batched_candidate_dict]
    except:
        batches_to_do = list(range(len(batched_ligand_pairs)))

    prefix = "Resuming" if batches_to_do[0] else "Starting"
    print(f"{prefix} job {job_name} from batch {batches_to_do[0]}/{len(batches_to_do)}...")

    batched_candidate_dict = defaultdict(dict)
    for curr_batch in progressbar(batches_to_do):
        with ProcessPoolExecutor(max_workers=num_cores) as executor:
            candidates = list(executor.map(get_candidates,
                                      batched_ligand_pairs[curr_batch]))
        """
                                      # Batch the batched list according to the
                                      # number of cores to use.
                                      get_batched_list(batched_ligand_pairs[curr_batch],
                                                       num_cores))
        candidates = [candidate_smiles for sublist in candidates
                      for candidate_smiles in sublist]
        """

        for i, smiles_pair in enumerate(batched_ligand_pairs[curr_batch]):
            batched_candidate_dict[curr_batch][smiles_pair]= candidates[i]

        with open(f"data/cache/{job_name}/batched_candidate_dict.pkl", "wb") as f:
            dump(batched_candidate_dict, f)

    candidate_dict = dict()
    for subdict in batched_candidate_dict.values():
        for smiles_pair, decoy_candidates in subdict.items():
            candidate_dict[smiles_pair] = decoy_candidates

    with open(f"data/candidate_dicts/{job_name}.pkl", "wb") as f:
        dump(candidate_dict, f)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Supply input and save paths.")
    parser.add_argument("job_name",
                        help=(f"Job name to be used to set the results cache "
                              f"location and decoy candidates."))
    # TODO: remove default after testing and make required=True.
    parser.add_argument("-i",
                        default="data/test_actives_protonated",
                        help=(f"Relative path to input actives_protonated "
                              f"SMILES file as created by 04_get_protonated.py"))
    parser.add_argument("--batch_size", default=128, type=int, help=BATCH_HELP_MESSAGE)
    parser.add_argument("--num_cores", type=int, default=multiprocessing.cpu_count(),
                        help=(f"Number of CPU cores to use for processing every"
                              f" batch Default: all available cores."))
    parser.add_argument("--reset", action="store_true", help=RESET_HELP_MESSAGE)
    args = vars(parser.parse_args())

    cache_dir = f"data/cache/{args['job_name']}"
    if not os.path.isdir(cache_dir):
        os.makedirs(cache_dir)
    elif args["reset"]:
        rmtree(cache_dir, ignore_errors=True)
        os.makedirs(cache_dir)

    print(f"Creating decoy candidates for {args['i']}")

    SORTED_ZINC_DFS = [create_sa(prop) for prop in COLUMNS[1:]]

    print("[1/2] Loading data...")
    batched_ligand_pairs = get_batched_ligand_pairs(args["i"], args["batch_size"])

    print(f"[2/2] Getting decoy candidates...")
    write_candidate_dict(batched_ligand_pairs,
                         args["num_cores"],
                         args["job_name"])

    print("Done.")

