"""Lorem ipsum."""

import os
import pickle
import argparse
import numpy as np
import pandas as pd
from rdkit.Chem import AllChem
from collections import Counter
from rdkit import Chem, DataStructs
from progressbar import progressbar
from collections import defaultdict


np.random.seed(12345)

zinc_smiles = np.load("data/all_zinc_smiles_ref.npy",
                      allow_pickle=True).flatten()

DATASET = "DUDE"


def get_processed_cand_dict(pdb_id):
    candidate_dict = pd.read_pickle(f"data/candidates/{pdb_id}_candidate_dict.pkl")
    # XXX: Consider if this can be done in 05.
    return  {(pdb_id, smiles): zinc_smiles[indices.astype(int)].tolist()
             for smiles, indices in candidate_dict.items()}


def get_fingerprint(smiles):
    return AllChem.GetMorganFingerprint(
                Chem.MolFromSmiles(smiles),
                2)


def get_sim(fprint_1, fprint_2):
    # TODO: Add choice of fingerprinting and (thus) similarity method.
    return DataStructs.DiceSimilarity(fprint_1, fprint_2)


def get_max_tc(cand_smiles, pairs):
    cand_fingerprint = get_fingerprint(cand_smiles)
    pairs_fingerprints = [get_fingerprint(pair[1]) for pair in pairs]
    return max([get_sim(cand_fingerprint, curr_fingerprint)
                for curr_fingerprint in pairs_fingerprints])


print(f"Writing decoys for {DATASET} dataset.")


print("[1/6] Concatenating decoy candidate dictionaries...")
all_pdb_ids = list(pd.read_pickle(
                    f"data/{DATASET.lower()}_actives_protonated_dict.pkl").keys())
all_candidate_dict = dict()
for pdb_id in progressbar(all_pdb_ids):
    all_candidate_dict.update(get_processed_cand_dict(pdb_id))

all_candidates = [smiles for sublist in all_candidate_dict.values()
                  for smiles in sublist]
candidate_counts = Counter(all_candidates)
raw_counts = list(candidate_counts.values())

print(f"{len(candidate_counts)} unique candidates out of a total of"
      f" {len(all_candidates)} "
      f"({len(candidate_counts) * 100 / len(all_candidates):.2f}%)")

print(f"Mean duplication: {np.mean(raw_counts):.0f} +- {np.std(raw_counts):.0f} "
      f"(min: {np.min(raw_counts):.0f}, max: {np.max(raw_counts):.0f})")


print("[2/6] Reversing concatenated dictionary...")
reversed_dict = defaultdict(list)
for pair in progressbar(all_candidate_dict):
    for cand_smiles in all_candidate_dict[pair]:
        reversed_dict[cand_smiles].append(pair)


print("[3/6] Getting maximum Tc values for each candidate...")
# TODO: if batched properly, this could be done many times faster.
if not os.path.exists("maximum_tc.pkl"):
    maximum_tc = defaultdict(int)
    repeat_fprints = dict()
    for pair in progressbar(all_candidate_dict):
        active_fingerprint = get_fingerprint(pair[1])
        for cand_smiles in all_candidate_dict[pair]:
            if candidate_counts[cand_smiles] - 1:
                try:
                    cand_fingerprint = repeat_fprints[cand_smiles]
                except:
                    cand_fingerprint = get_fingerprint(cand_smiles)
                    repeat_fprints[cand_smiles] = cand_fingerprint

                maximum_tc[cand_smiles] = max(maximum_tc[cand_smiles],
                                              DataStructs.DiceSimilarity(
                                                    active_fingerprint,
                                                    cand_fingerprint))
            else:
                cand_fingerprint = get_fingerprint(cand_smiles)
                maximum_tc[cand_smiles] = DataStructs.DiceSimilarity(
                                                    active_fingerprint,
                                                    cand_fingerprint)
    with open("maximum_tc.pkl", "wb") as f:
        pickle.dump(maximum_tc, f)
else:
    maximum_tc = pd.read_pickle("maximum_tc.pkl")


print("[4/6] Getting qualifying candidate decoys...")
# XXX: What's the purpose of this relative cutoff relative to a hard one?
maximum_tc = [[cand_smiles, max_tcs]
              for cand_smiles, max_tcs in maximum_tc.items()]
maximum_tc = sorted(maximum_tc, key=lambda pair: pair[1])
maximum_tc = maximum_tc[:(len(maximum_tc) // 4)]
qualifying_cands = dict(maximum_tc)
max_tcs = list(qualifying_cands.values())
print(f"Got {len(qualifying_cands)} for {len(all_candidate_dict)} actives.")
print(f"Mean Tc: {np.mean(max_tcs):.3f} +- {np.std(max_tcs):.3f} "
      f"(min: {np.min(max_tcs):.3f}, max: {np.max(max_tcs):.3f})")

all_candidate_dict = {active_smiles: [cand for cand in sublist
                                       if cand in qualifying_cands]
                       for active_smiles, sublist in all_candidate_dict.items()}
num_active_smiles = len(all_candidate_dict)
all_candidate_dict = {active_smiles: sublist for active_smiles, sublist in
                       all_candidate_dict.items() if len(sublist) >= 50}
print(f"{100 - len(all_candidate_dict) * 100 / num_active_smiles:.2f}% actives "
      f"do not have at least 50 qualifying candidate decoys; not assigning any to them.")

print("[5/6] Assigning qualifying candidates...")
for pair, sublist in progressbar(all_candidate_dict.items()):
    all_candidate_dict[pair] = \
        sorted(sublist,
               key=lambda cand_smiles: candidate_counts[cand_smiles])[:50]

# Print something about duplication here.
all_decoys = [decoy_smiles for sublist in all_candidate_dict.values()
              for decoy_smiles in sublist]
print(f"{len(set(all_decoys)) * 100 / len(all_decoys):.2f}% of assigned decoys are unique.")

print("[6/6] Writing assigned decoys to disk...")
rows = []
for pair, sublist in progressbar(all_candidate_dict.items()):
    for decoy_smiles in sublist:
        rows.append(f"{pair[0]} {pair[1]} {decoy_smiles}")
np.random.shuffle(rows)

with open(f"data/{DATASET.lower()}_decoys", "w") as f:
    f.write("\n".join(rows))

print("Done.")

