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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Supply job name.")
    parser.add_argument("job_name",
                        help="Job name to be used to load the decoy candidates.")
    args = vars(parser.parse_args())

    print("Selecting decoys from qualifying candidates for each active...")

    print(f"[1/X] Loading decoy candidate dictionary...")
    candidate_dict = \
        {smiles_pair: ZINC_SMILES[indices.astype(int)].tolist()
         for smiles_pair, indices in
         pd.read_pickle(f"data/candidate_dicts/{args['job_name']}.pkl").items()
         if len(indices)}

    with open(f"data/cache/{args['job_name']}/qualifying_candidates", "r") as f:
        qualifying_candidates = [smiles.strip("\n") for smiles in f.readlines()]
    qualifying_candidates = dict(zip(qualifying_candidates, [[] for _ in
                                                             qualifying_candidates]))

    for smiles_pair in progressbar(candidate_dict):
        candidate_dict[smiles_pair] = [curr_smiles for curr_smiles in
                                    candidate_dict[smiles_pair] if curr_smiles in
                                       qualifying_candidates]
    candidate_dict = {smiles_pair: candidate_smiles for smiles_pair,
candidate_smiles in candidate_dict.items()
                      if len(candidate_smiles)}

    all_qualified = [smiles for candidate_smiles in candidate_dict.values()
                     for smiles in candidate_smiles]
    candidate_counts = Counter(all_qualified)

    sorting_key = lambda x: candidate_counts[x]

    for smiles_pair, candidate_smiles in progressbar(candidate_dict.items()):
        if len(candidate_smiles) < 50:
            continue
        candidate_dict[smiles_pair] = sorted(candidate_smiles,
                                            key= sorting_key)[:50]

    actives_dict = defaultdict(set)
    failed_actives = []
    for smiles_pair, candidate_smiles in progressbar(candidate_dict.items()):
        actives_dict[smiles_pair.split("_")[0]].update(candidate_smiles)

    for active_smiles, candidate_smiles in progressbar(actives_dict.items()):
        if len(candidate_smiles) < 50:
            failed_actives.append(active_smiles)
        else:
            actives_dict[active_smiles] = \
                sorted(list(candidate_smiles), key=sorting_key)[:50]

    actives_dict = {active_smiles: candidate_smiles for active_smiles,
candidate_smiles in actives_dict.items() if len(candidate_smiles) >= 50}

    print(len(failed_actives), "failed")

    to_write = []
    for active_smiles, decoy_list in actives_dict.items():
        for decoy_smiles in decoy_list:
            to_write.append(f"{active_smiles} {decoy_smiles}")

    with open(f"data/{args['job_name']}_decoys", "w") as f:
        f.write("\n".join(to_write))

    print("Done.")

