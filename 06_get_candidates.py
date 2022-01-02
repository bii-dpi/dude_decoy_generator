"""Find 3000-9000 decoy candidates for each active (and its protonated states)."""

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

print("Finding decoy candidates for each active (and its protonated states)...")


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


def get_smiles_in_range(value, delta, curr_df):
    """Get the SMILES indices that are within the property value range."""

    # XXX: This can perhaps be made very slightly tighter.
    left_pos = np.searchsorted(curr_df[1], value - delta, side="left")
    right_pos = np.searchsorted(curr_df[1], value + delta, side="right")

    if left_pos == len(curr_df[1]):
        left_pos -= 1
    if right_pos == len(curr_df[1]):
        right_pos -= 1

    assert right_pos >= left_pos

    return curr_df[0][left_pos: right_pos + 1]


def get_property_matched(curr_ligand_props, window_index):
    """Get the SMILES indices that are in-range for all six properties."""

    # Keep the candidates that fall within property ranges for all six
    # properties.
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


def get_candidates(smiles_pair):
    """Get candidate decoy indices for a single ligand."""

    smiles = smiles_pair.split("_")[1]
    mol = Chem.MolFromSmiles(smiles)
    ligand_props = get_properties(mol)

    prop_matched_indices = []
    window_index = 0
    while len(prop_matched_indices) < 9000 and window_index < 7:
        # Get the decoy candidate indices that correspond to the current range.
        latest = get_property_matched(ligand_props,
                                      window_index)
        window_index += 1

        # If there were no candidates available for this strict a range,
        # continue to the next iteration with a wider range.
        if latest is None:
            continue

        # If there are no candidate indices currently stored, we can potentially
        # use all of the latest. Otherwise, we will take a union of the two.
        if not len(prop_matched_indices):
            tried_union = latest.copy()
        else:
            tried_union = np.union1d(prop_matched_indices, latest)

        # If this attempted new indices set is too large, it must be reduced in
        # size. How exactly this is done depends on the number of novel indices
        # that were found in the current iteration.
        if len(tried_union) > 9000:
            latest = np.setdiff1d(latest, prop_matched_indices,
                                  assume_unique=True)
            # We try to sample exactly the number of new indices we need from
            # the latest found to create a set of 9000. If this does not work
            # because we have too few newly-found, we can include all of the
            # novel indices into the indices set.
            try:
                prop_matched_indices = \
                    np.union1d(prop_matched_indices,
                               np.random.choice(latest,
                                                9000 - len(prop_matched_indices),
                                                replace=False))
            except Exception as e:
                prop_matched_indices = np.union1d(latest,
                                                  prop_matched_indices)
        else:
            prop_matched_indices = tried_union.copy()

        # If we have at least 3000 suitable indices, we do not need to try to
        # find more.
        if len(prop_matched_indices) >= 3000:
            break

    return prop_matched_indices


def write_candidate_dict(batched_ligand_pairs, num_cores, job_name):
    """Write the candidate dictionary for the dataset."""

    # Find the indices of the batches still to be completed.
    try:
        batched_candidate_dict = \
            pd.read_pickle(f"data/cache/{job_name}/batched_candidate_dict.pkl")
        batches_to_do = [i for i in range(len(batched_ligand_pairs))
                         if i not in batched_candidate_dict]
    except:
        batched_candidate_dict = defaultdict(dict)
        batches_to_do = list(range(len(batched_ligand_pairs)))

    prefix = "Resuming" if batches_to_do[0] else "Starting"
    print(f"{prefix} job {job_name} from batch {batches_to_do[0]}/{len(batches_to_do)}")

    for curr_batch in progressbar(batches_to_do):
        # Get the list of candidates for each SMILES pair in the batch.
        with ProcessPoolExecutor(max_workers=num_cores) as executor:
            candidates = list(executor.map(get_candidates,
                                           batched_ligand_pairs[curr_batch]))

        # Store these lists in the batched candidte dictionary, under this batch
        # number.
        for i, smiles_pair in enumerate(batched_ligand_pairs[curr_batch]):
            batched_candidate_dict[curr_batch][smiles_pair] = candidates[i]

        # Overwrite the batched candidate dictionary.
        with open(f"data/cache/{job_name}/batched_candidate_dict.pkl", "wb") as f:
            dump(batched_candidate_dict, f)

    # Now that all batches have been processed, remove the batch number from the
    # batched candidate dictionary to create the candidate dictionary.
    candidate_dict = dict()
    for subdict in batched_candidate_dict.values():
        candidate_dict.update(subdict)

    with open(f"data/candidate_dicts/{job_name}.pkl", "wb") as f:
        dump(candidate_dict, f)


def get_batched_ligand_pairs(input_path, batch_size):
    """Get referenced ligand pairs list split into batches."""

    with open(input_path, "r") as f:
       ligand_pairs = [line.strip("\n") for line in f.readlines()]

    batched_ligand_pairs = []
    indices = list(range(0, len(ligand_pairs), batch_size)) + [-1]
    for i in range(len(indices) - 1):
        if indices[i + 1] == -1:
            batched_ligand_pairs.append(ligand_pairs[indices[i]:])
        else:
            batched_ligand_pairs.append(ligand_pairs[indices[i]:indices[i + 1]])

    return batched_ligand_pairs


def create_sa(prop):
    """Put the pre-sorted SMILES arrays and properties in shared memory."""

    smiles = np.load(f"data/sorted/{prop}_smiles_reorder.npy")
    values = np.load(f"data/sorted/{prop}_values.npy")

    smiles_ = sa.create(f"file://data/sa/{prop}_smiles", smiles.shape, dtype=smiles.dtype)
    values_ = sa.create(f"file://data/sa/{prop}_values", values.shape, dtype=values.dtype)

    np.copyto(smiles_, smiles)
    np.copyto(values_, values)

    return smiles_, values_


if __name__ == "__main__":
    SORTED_ZINC_DFS = [create_sa(prop) for prop in COLUMNS[1:]]

    parser = argparse.ArgumentParser(description="Supply job name and input path.")
    parser.add_argument("job_name",
                        help=("Job name to be used to set the results cache "
                              "location and save the decoy candidates."))
    # TODO: remove default after testing and make required=True.
    parser.add_argument("-i",
                        default="data/test_actives_protonated",
                        help=("Relative path to input actives_protonated "
                              "SMILES file as created by 04_get_protonated.py"))
    parser.add_argument("--batch_size", default=128, type=int, help=BATCH_HELP_MESSAGE)
    parser.add_argument("--num_cores", type=int, default=multiprocessing.cpu_count(),
                        help=("Number of CPU cores to use for processing every"
                              " batch Default: all available cores."))
    parser.add_argument("--reset", action="store_true", help=RESET_HELP_MESSAGE)
    args = vars(parser.parse_args())

    cache_dir = f"data/cache/{args['job_name']}"
    if not os.path.isdir(cache_dir):
        os.makedirs(cache_dir)
    elif args["reset"]:
        rmtree(cache_dir, ignore_errors=True)
        os.makedirs(cache_dir)

    print("[1/2] Loading data...")
    batched_ligand_pairs = get_batched_ligand_pairs(args["i"], args["batch_size"])

    print(f"[2/2] Getting decoy candidates...")
    write_candidate_dict(batched_ligand_pairs,
                         args["num_cores"],
                         args["job_name"])

    print("Done.")


