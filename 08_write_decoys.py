"""Write decoys from qualifying candidates to disk."""

import argparse
import numpy as np
import pandas as pd
from progressbar import progressbar
from collections import Counter, defaultdict


print("Selecting and writing decoys for each active...")


ZINC_SMILES = np.load("data/all_zinc_smiles_ref.npy",
                      allow_pickle=True).flatten()
ID_dict = pd.read_pickle("data/zinc_dict.pkl")


def write_decoys(path, dict_):
    """Write decoys in the input dictionary to path."""

    to_write = []
    for active_smiles, decoy_list in dict_.items():
        for decoy_smiles in decoy_list:
            to_write.append(f"{active_smiles} {decoy_smiles} {ID_dict[decoy_smiles]}")

    with open(path, "w") as f:
        f.write("\n".join(to_write))


def make_percentage(freq, counts_sum):
    """Get percentage as a pretty string."""

    perc = (freq * 100) // counts_sum
    if perc < 1:
        perc = "<1"

    return f"{perc}%"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Supply job name and output path.")
    parser.add_argument("job_name",
                        help="Job name to be used to load the decoy candidates.")
    parser.add_argument("-o", default=None,
                        help=("Decoy file output path; do not end it with an extension."
                              " Default: data/\{job_name\}_decoys"))
    args = vars(parser.parse_args())

    output_path = args["o"]
    if not output_path:
        output_path = f"data/{args['job_name']}_decoys"
        print(f"No output path supplied; writing decoys to {output_path}")

    print(f"[1/3] Loading decoy candidate dictionary...")
    candidate_dict = \
        {smiles_pair: ZINC_SMILES[indices.astype(int)].tolist()
         for smiles_pair, indices in
         pd.read_pickle(f"data/candidate_dicts/{args['job_name']}.pkl").items()
         if len(indices)}

    # Load qualifying candidates, as generated by 07, to filter candidate_dict.
    with open(f"data/cache/{args['job_name']}/qualifying_candidates", "r") as f:
        qualifying_candidates = [smiles.strip("\n") for smiles in f.readlines()]
    # Create a hash table for quick look-up.
    qualifying_candidates = dict(zip(qualifying_candidates,
                                     [[] for _ in qualifying_candidates]))

    print("Filtering candidates")
    for smiles_pair in progressbar(candidate_dict):
        candidate_dict[smiles_pair] = \
            [curr_smiles for curr_smiles in candidate_dict[smiles_pair]
             if curr_smiles in qualifying_candidates]

    # Remove any SMILES pairs without any candidates left after filtration.
    candidate_dict = {smiles_pair: candidate_smiles
                      for smiles_pair, candidate_smiles in candidate_dict.items()
                      if len(candidate_smiles)}

    # Consolidate candidates across the protonated states of each active.
    print("Consolidating candidates by active")
    decoys_dict = defaultdict(set)
    for smiles_pair, candidate_smiles in progressbar(candidate_dict.items()):
        decoys_dict[smiles_pair.split("_")[0]].update(candidate_smiles)

    print(f"[2/3] Assigning decoys to actives by duplicatedness...")
    # Calculate candidate counts.
    all_qualified = [smiles for candidate_smiles in candidate_dict.values()
                     for smiles in candidate_smiles]
    candidate_counts = Counter(all_qualified)

    failed_actives = []
    for active_smiles, candidate_smiles_set in progressbar(decoys_dict.items()):
        if len(candidate_smiles_set) < 50:
            failed_actives.append(active_smiles)

        decoys_dict[active_smiles] = \
            sorted(list(candidate_smiles_set),
                   key=lambda smiles: candidate_counts[smiles])[:50]

    if len(failed_actives):
        failed_output_path = f"{output_path}_failed"
        print(f"{len(failed_actives)}/{len(decoys_dict)} "
              f"({make_percentage(len(failed_actives), len(decoys_dict))}) "
              f"actives have fewer than 50 assignable decoys; writing their "
              f"decoys to a separate file: {failed_output_path}")
        failed_dict = dict()
        for failed_active in failed_actives:
            failed_dict[failed_active] = decoys_dict[failed_active]
            del decoys_dict[failed_active]

    # Calculate assigned decoy statistics.
    all_decoys = [smiles for decoy_smiles_list in decoys_dict.values()
                  for smiles in decoy_smiles_list]
    decoy_counts = Counter(list(Counter(all_decoys).values()))
    counts_sum = np.sum(list(decoy_counts.values()))
    decoy_props = {count: make_percentage(freq, counts_sum)
                   for count, freq in decoy_counts.items()}
    pretty_props = ""
    for count in sorted(list(decoy_props.keys())):
        pretty_props += f"Assigned to {count} actives: {decoy_props[count]}\n"

    print(f"{len(set(all_decoys))} unique decoys assigned")
    print(pretty_props[:-1])

    print(f"[3/3] Writing decoys...")
    write_decoys(output_path, decoys_dict)
    write_decoys(failed_output_path, failed_dict)

    print("Done.")

