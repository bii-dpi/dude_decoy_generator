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
import pickle
import argparse
import numpy as np
import pandas as pd
from rdkit.Chem import AllChem
from rdkit import Chem, DataStructs
from progressbar import progressbar
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor


np.random.seed(12345)


ZINC_SMILES = np.load("data/all_zinc_smiles_ref.npy",
                      allow_pickle=True).flatten()


def get_fingerprint(smiles):
    """Get the ECFP4 fingerprint of the input SMILES' compound."""

    return AllChem.GetMorganFingerprint(
                Chem.MolFromSmiles(smiles),
                2)


def get_fprint_subdict(input_subset, i):
    """Lorem ipsum."""

    return dict(zip(input_subset,
                    pd.read_pickle(f"data/fingerprints/ecfp4/{i}.pkl")))


def get_sim(fprint_1, fprint_2):
    """Lorem ipsum."""

    # TODO: Add choice of fingerprinting and (thus) similarity method.
    return DataStructs.DiceSimilarity(fprint_1, fprint_2)


def get_maximum_tc(triplet):
    """Lorem ipsum."""

    smiles, cand_fprint, active_fprints = triplet

    return smiles, max([get_sim(cand_fprint, active_fprint)
                        for active_fprint in active_fprints])


def get_maximum_tc_batched(triplets):
    """Lorem ipsum."""

    return [get_maximum_tc(triplet) for triplet in triplets]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Supply job name.")
    parser.add_argument("job_name",
                        help="Job name to be used to load the decoy candidates.")
    args = vars(parser.parse_args())

    print("Selecting decoys from candidates for each active...")

    print(f"[1/X] Loading decoy candidate dictionary...")
    candidate_dict = \
        {smiles_pair: ZINC_SMILES[indices.astype(int)].tolist()
         for smiles_pair, indices in
         pd.read_pickle(f"data/candidate_dicts/{args['job_name']}.pkl").items()
         if len(indices)}

    print(f"[2/X] Getting active and decoy candidate fingerprints...")

    # Get actives' fingerprints.
    print("Actives")
    actives_fprint_dict = dict()
    for smiles_pair in progressbar(candidate_dict):
        actives_fprint_dict[smiles_pair] = \
            get_fingerprint(smiles_pair.split("_")[1])

    with open("data/all_zinc.smi", "r") as f:
        all_zinc = [smiles.strip("\n") for smiles in f.readlines()]

    # Get decoy candidates' fingerprints.
    print("Decoy candidates")
    candidates_fprint_dict = dict()
    indices = list(range(0, len(all_zinc), 10000)) + [-1]
    input_subsets = []
    for i in range(len(indices) - 1):
        if indices[i + 1] != -1:
            input_subsets.append((all_zinc[indices[i]:indices[i + 1]], indices[i]))
        else:
            input_subsets.append((all_zinc[indices[i]:], indices[i]))

    for input_subset, i in progressbar(input_subsets):
        candidates_fprint_dict.update(get_fprint_subdict(input_subset, i))

    print(f"[3/X] Calculating maximum Tc between decoy candidates and "
          f"corresponding actives...")
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

    with ProcessPoolExecutor() as executor:
        maximum_tc = executor.map(get_maximum_tc_batched, reversed_list_batched)
    maximum_tc = [pair for sublist in maximum_tc for pair in sublist]

    print("[4/X] Sorting candidate SMILES by maximum Tc to keep the X% with the lowest Tcs...")
    # Should make this variable.
    maximum_tc = sorted(maximum_tc, key=lambda pair: pair[1])
    with open(f"data/cache/{args['job_name']}/qualifying_candidates", "w") as f:
        f.write("\n".join(maximum_tc[:int(len(maximum_tc) / 4)]))

    print("Done.")

