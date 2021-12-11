import os
import pickle

import numpy as np
import pandas as pd

from active import Active
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

print("[1/X] Concatenating decoy candidate dictionaries...")
all_pdb_ids = list(pd.read_pickle(
                    f"data/{DATASET.lower()}_actives_protonated_dict.pkl").keys())
all_candidate_dict = dict()
for pdb_id in progressbar(all_pdb_ids):
    all_candidate_dict.update(get_processed_cand_dict(pdb_id))

all_candidates = [smiles for sublist in all_candidate_dict.values()
                  for smiles in sublist]
candidate_counts = Counter(all_candidates)

print(len(candidate_counts) / len(all_candidates))


print("[2/X] Reversing concatenated dictionary...")
reversed_dict = defaultdict(list)
for pair in progressbar(all_candidate_dict):
    for cand_smiles in all_candidate_dict[pair]:
        reversed_dict[cand_smiles].append(pair)


print("[3/X] Getting maximum Tc values for each candidate...")
"""
reversed_dict = {cand_smiles: sublist for cand_smiles, sublist in
reversed_dict.items() if cand_smiles in qualifying_cands}

print(len(reversed_dict) / 102)

# XXX: Much slower to do it in this more intuitive way.
maximum_tc = dict()
freqs = dict()
for cand_smiles, pairs in progressbar(reversed_dict.items()):
    maximum_tc[cand_smiles] = get_max_tc(cand_smiles, pairs)
    freqs[cand_smiles] = len(pairs)


"""
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


print("[4/X] Getting qualifying candidate decoys...")
maximum_tc = [[cand_smiles, max_tcs]
              for cand_smiles, max_tcs in maximum_tc.items()]
maximum_tc = sorted(maximum_tc, key=lambda pair: pair[1])
maximum_tc = maximum_tc[:(len(maximum_tc) // 4)]
qualifying_cands = dict(maximum_tc)
max_tcs = list(qualifying_cands.values())
print(f"Got {len(qualifying_cands)} for {len(reversed_dict)} actives.")
print(f"Mean Tc: {np.mean(max_tcs):.3f} +- {np.std(max_tcs):.3f} "
      f"(min: {np.min(max_tcs):.3f}, max: {np.max(max_tcs):.3f})")

all_candidate_dict = {active_smiles: [cand for cand in sublist
                                       if cand in qualifying_cands]
                       for active_smiles, sublist in all_candidate_dict.items()}
num_active_smiles = len(all_candidate_dict)
all_candidate_dict = {active_smiles: sublist for active_smiles, sublist in
                       all_candidate_dict.items() if len(sublist) >= 50}
print(f"{len(all_candidate_dict) * 100 / num_active_smiles:.2f}% actives do not have "
      f"at least 50 qualifying candidate decoys; not assigning any to them.")

print("[5/X] Assigning qualifying candidates...")
actives = [Active(active_pair, assignable_decoys)
                   for active_pair, assignable_decoys in
                   all_candidate_dict.items()]

reversed_dict = {cand_smiles: actives_smiles
                 for cand_smiles, actives_smiles in reversed_dict.items()
                 if cand_smiles in qualifying_cands}

qualifying_cands = [[cand_smiles, len(actives_smiles)]
                    for cand_smiles, actives_smiles in reversed_dict.items()]
qualifying_cands = [pair[0] for pair in
                    sorted(qualifying_cands, key=lambda pair: pair[1])]


num_done = 0
while not num_done == len(all_candidate_dict):
    for cand in qualifying_cands[:3000]:
        for active in actives:
            if not active.has_50_decoys() and active.can_assign(cand):
                active.add_decoy(cand)
                break
        actives.sort()

    num_done = np.sum([active.has_50_decoys() for active in actives])
    print(f"{num_done * 100 / len(all_candidate_dict):.2f}% done.",
          end="\r")
print("")


print("[6/X] Writing assigned decoys to disk...")

print("Done. (and some stats)")
