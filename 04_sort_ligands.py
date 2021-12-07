import numpy as np

import pandas as pd
from os import remove
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor


COLUMNS = ("smiles", "mol_wt", "logp",
           "num_rotatable", "num_hba", "num_hbd", "net_charge")

all_zinc_df = pd.read_csv("data/all_zinc_w_prop.csv",
                          names=COLUMNS,
                          dtype={"smiles": str, "mol_wt": np.float32,
                                 "logp": np.float32, "num_rotatable": int,
                                 "num_hba": int, "num_hbd": int, "net_charge": int})
print(all_zinc_df.columns)

for column in tqdm(COLUMNS[1:]):
    all_zinc_df[["smiles", column]].\
        sort_values(by=column).\
        to_csv(f"data/all_zinc_w_prop_{column}.csv", index=False)

