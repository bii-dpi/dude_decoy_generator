import pickle

import numpy as np
import pandas as pd

from rdkit.Chem import AllChem
from collections import Counter
from rdkit import Chem, DataStructs
from progressbar import progressbar
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor


DATASET = "DUDE"


zinc_smiles = np.load("data/all_zinc_smiles_ref.npy",
                      allow_pickle=True).flatten()


def get_processed_cand_dict(pdb_id):
    candidate_dict = pd.read_pickle(f"data/candidates/{pdb_id}_candidate_dict.pkl")
    # XXX: Consider if this can be done in 05.
    return  {(pdb_id, smiles): zinc_smiles[indices.astype(int)].tolist()
             for smiles, indices in candidate_dict.items()}


def get_fingerprint(smiles):
    return AllChem.GetMorganFingerprint(
                Chem.MolFromSmiles(smiles),
                2)


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
for pair in all_candidate_dict:
    for cand_smiles in all_candidate_dict[pair]:
        reversed_dict[cand_smiles].append(pair)


print("[3/X] Getting maximum Tc values for each candidate...")
# TODO: if batched properly, this could be done many times faster.
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

maximum_tc = [[cand_smiles, max_tcs]
              for cand_smiles, max_tcs in maximum_tc.items()]
maximum_tc = sorted(maximum_tc, key=lambda pair: pair[1])
maximum_tc = maximum_tc[:(len(maximum_tc) // 4)]


