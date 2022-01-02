"""Filter candidate decoys by maximum Tc."""

import os
import pickle
import argparse
import numpy as np
import pandas as pd
from rdkit.Chem import AllChem
from rdkit import Chem, DataStructs
from progressbar import progressbar
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor


print("Filtering candidates by Tc...")

np.random.seed(12345)

ZINC_SMILES = np.load("data/all_zinc_smiles_ref.npy",
                      allow_pickle=True).flatten()


def get_fingerprint(smiles):
    """Get the ECFP4 fingerprint of the input SMILES' compound."""

    # There is no need to check if the input SMILES can be converted to a Mol
    # successfully, as previous processing steps have removed SMILES that cannot
    # be.
    return AllChem.GetMorganFingerprint(Chem.MolFromSmiles(smiles), 2)


def get_fprint_subdict(input_subset, i):
    """Load the ith ECFP4 fingerprint dictionary fragment."""

    return dict(zip(input_subset,
                    pd.read_pickle(f"data/fingerprints/ecfp4/{i}.pkl")))


def get_sim(fprint_1, fprint_2):
    """Get the Tc between two ECFP4 fingerprints."""

    return DataStructs.DiceSimilarity(fprint_1, fprint_2)


def get_maximum_tc(triplet):
    """Get the maximum Tc between the decoy candidate and its actives."""

    cand_smiles, cand_fprint, active_fprints = triplet

    return cand_smiles, max([get_sim(cand_fprint, active_fprint)
                             for active_fprint in active_fprints])


def get_maximum_tc_batched(triplets):
    """Get the maximum Tc associated with each decoy candidate in the batch."""

    return [get_maximum_tc(triplet) for triplet in triplets]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Supply job name.")
    parser.add_argument("job_name",
                        help="Job name to be used to load the decoy candidates.")
    args = vars(parser.parse_args())

    print(f"[1/4] Loading decoy candidate dictionary...")
    candidate_dict = \
        {smiles_pair: ZINC_SMILES[indices.astype(int)].tolist()
         for smiles_pair, indices in
         pd.read_pickle(f"data/candidate_dicts/{args['job_name']}.pkl").items()
         if len(indices)}

    print(f"[2/4] Getting active and decoy candidate fingerprints...")
    # Get actives' fingerprints.
    print("Actives")
    actives_fprint_dict = dict()
    for smiles_pair in progressbar(candidate_dict):
        actives_fprint_dict[smiles_pair] = \
            get_fingerprint(smiles_pair.split("_")[1])

    # Get decoy candidates' fingerprints.
    print("Decoy candidates")
    with open("data/all_zinc.smi", "r") as f:
        all_zinc = [smiles.strip("\n") for smiles in f.readlines()]

    indices = list(range(0, len(all_zinc), 10000)) + [-1]
    input_subsets = []
    for i in range(len(indices) - 1):
        if indices[i + 1] != -1:
            input_subsets.append((all_zinc[indices[i]:indices[i + 1]], indices[i]))
        else:
            input_subsets.append((all_zinc[indices[i]:], indices[i]))

    candidates_fprint_dict = dict()
    for input_subset, i in progressbar(input_subsets):
        candidates_fprint_dict.update(get_fprint_subdict(input_subset, i))

    print(f"[3/4] Calculating maximum Tc between decoy candidates and "
          f"corresponding actives...")
    print("Preparing data")
    # Reverse the candidate dictionary to allow parallelization of maximum Tc
    # calculation.
    reversed_dict = defaultdict(list)
    for smiles_pair in progressbar(candidate_dict):
        for cand_smiles in candidate_dict[smiles_pair]:
            reversed_dict[cand_smiles].append(smiles_pair)

    reversed_list = [[cand_smiles, candidates_fprint_dict[cand_smiles],
                     [actives_fprint_dict[actives_smiles_pair]
                      for actives_smiles_pair in actives_smiles_pair_list]]
                     for cand_smiles, actives_smiles_pair_list
                     in reversed_dict.items()]

    indices = list(range(0, len(reversed_list), 10000)) + [-1]
    reversed_list_batched = []
    for i in range(len(indices) - 1):
        if indices[i + 1] == -1:
            reversed_list_batched.append(reversed_list[indices[i]:])
        else:
            reversed_list_batched.append(reversed_list[indices[i]:indices[i + 1]])

    print("Calculating")
    with ProcessPoolExecutor() as executor:
        maximum_tc = executor.map(get_maximum_tc_batched, reversed_list_batched)
    maximum_tc = [pair for sublist in maximum_tc for pair in sublist]

    print("[4/4] Sorting candidates by maximum Tc to keep the bottom 25%...")
    maximum_tc = sorted(maximum_tc, key=lambda pair: pair[1])
    maximum_tc = [pair[0] for pair in maximum_tc[:(len(maximum_tc) // 4)]]
    with open(f"data/cache/{args['job_name']}/qualifying_candidates", "w") as f:
        f.write("\n".join(maximum_tc))

    print("Done.")

