"""Generate Tc matrix."""


import os
import pickle

import numpy as np
import pandas as pd

from rdkit import Chem
from itertools import combinations
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor


SAVE_MOLS = True
SAVE_FINGERPRINTS = True

SMILES_LIST_PATH = "example_smiles"
BATCH_SIZE = 10000


if __name__ == "__main__":
    print("Getting Mols...")
    if SAVE_MOLS:
        # Load SMILES list.
        with open(SMILES_LIST_PATH, "r") as f:
            all_smiles = f.readlines()

        # Precompute the Daylight fingerprints of each compound.
        # This is worth doing in parallel only if there are many
        # (tens of thousands) of SMILES: otherwise, the overhead
        # will make it slower to parallelize than run serially.

        # We batch the SMILES list into batches of 10000 to obtain
        # a good parallelizing overhead-speed tradeoff. Smaller batches
        # are also justifiable.
        indices = list(range(0, len(all_smiles), BATCH_SIZE)) + [-1]
        smiles_subsets = []
        for i in range(len(indices) - 1):
            if indices[i + 1] != -1:
                smiles_subsets.append(all_smiles[indices[i]: indices[i + 1]])
            else:
                smiles_subsets.append(all_smiles[indices[i]:])

        print("lala")
        # Create the Mol object list.
        with ProcessPoolExecutor() as executor:
            all_mols = executor.map(lambda smiles_subset: [Chem.MolFromSmiles(smiles)
                                                           for smiles in smiles_subset],
                                    smiles_subsets)

        # Flatten Mol object list.
        all_mols = [output for sublist in all_mols
                    for output in sublist]

        with open(f"{SMILES_LIST_PATH}_mols.pkl", "wb") as f:
            pickle.dump(all_mols, f)
    else:
        all_mols = pd.read_pickle(f"{SMILES_LIST_PATH}_mols.pkl")


    print("Getting fingerprints...")
    if SAVE_FINGERPRINTS:
        # As before, batch the Mols objects into subsets.
        indices = list(range(0, len(all_mols), BATCH_SIZE)) + [-1]
        mols_subsets = []
        for i in range(len(indices) - 1):
            if indices[i + 1] != -1:
                mols_subsets.append(all_fingerprints[indices[i]: indices[i + 1]])
            else:
                mols_subsets.append(all_fingerprints[indices[i]:])

        # Create the Fingerprints list.
        with ProcessPoolExecutor() as executor:
            all_fingerprints = executor.map(lambda mols_subset: [Chem.RDKFingerprint(smiles)
                                                                 for mol in mols_subset],
                                            mols_subsets)
        all_fingerprints = list(all_fingerprints)

        # Given the size of Daylight fingerprints as compared to Mol objects,
        # it's much faster to save and load this list in batches, without
        # flattening it.

        # Attach the batch index suffix to each subset.
        for i in range(len(indices) - 1):
            all_fingerprints[i] = (all_fingerprints[i], indices[i])

        def save_subset(input_tuple):
            with open(f"fingerprints/fingerprints_{input_tuple[1]}.pkl", "wb") as f:
                pickle.dump(input_tuple[0], f)

        # ThreadPoolExecutor() is better for IO-bound tasks.
        with ThreadPoolExecutor() as executor:
            executor.map(save_subset, all_fingerprints)

        all_fingerprints = [fingerprint for sublist in all_fingerprints
                            for fingerprint in sublist]
    else:
        fingerprint_fnames = [fname for fname in os.listdir("fingerprints/")
                              if fname.startswith("fingerprints_")]
        suffixes = [int(fname.split("_")[1].strip(".pkl"))
                    for fname in fingerprint_fnames]
        fingerprints = []
        # Loading fingerprints serially is quick enough; only saving them
        # requires parallelization.
        for suffix in sorted(suffixes):
            fingerprints += pd.read_pickle(f"fingerprints/fingerprints_{suffix}.pkl")

    print("Getting Tc matrix...")
    # Get Tc for all pairs of fingerprints.
    def get_tc(index_1):
        return {index_2: DataStructs.TanimotoSimilarity(fingerprints[index_1],
                                                        fingerprints[index_2])
                for index_2 in all_index_pairs[index_1]}


    index_combos = combinations(range(len(fingerprints)), 2)
    all_index_pairs = defaultdict(list)
    for index_1, index_2 in index_combos:
        all_index_pairs[index_1].append(index_2)

    with ProcessPoolExecutor() as executor:
        tc_matrix = executor.map(get_tc, range(len(fingerprints)))

    with open("tc_matrix.pkl", "wb") as f:
        pickle.dump(tc_matrix, f)


    print("Done.")

