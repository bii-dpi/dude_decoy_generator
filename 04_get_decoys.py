import numpy as np

from rdkit import Chem, DataStructs, AllChem


ATD_RATIO = 50


def get_sim(active_fingerprint, candidate_smiles):
    """Get Tc between the compounds' ECP4 fingerprints."""
    candidate_fingerprint = AllChem.GetMorganFingerprint(
                                Chem.MolFromSmiles(candidate_smiles),
                                2)
    return DataStructs.DiceSimilarity(active_fingerprint,
                                      candidate_fingerprint)


def get_curr_candidates(all_candidates, active_fingerprint):
    """Get the bottom 25% of candidates by Tc to the active."""
    all_sims = [get_sim(active_fingerprint, candidate_smiles)
                for candidate_smiles in all_candidates]
    return all_candidates[np.argsort(all_sims)][:int(len(all_candidates) * 0.25)]


def get_decoys(all_candidates, active_smiles):
    """Get 50 decoys for single SMILES."""
    print("Getting dissimilar-enough candidates...")
    active_fingerprint = AllChem.GetMorganFingerprint(
                            Chem.MolFromSmiles(active_smiles),
                            2)
    curr_candidates = get_curr_candidates(all_candidates, active_fingerprint)

    assert len(curr_candidates) >= 50, \
            f"{active_smiles} doesn't have enough dissimilar candidates."

