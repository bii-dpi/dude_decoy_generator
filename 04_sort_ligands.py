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
                                 "logp": np.float32, "num_rotatable": np.int16,
                                 "num_hba": np.int16, "num_hbd": np.int16,
                                 "net_charge": np.int16})
print(all_zinc_df.columns)

np.save("data/all_zinc_smiles_ref.npy", all_zinc_df[["smiles"]].to_numpy())

for column in tqdm(COLUMNS[1:]):
    smiles_reorder = all_zinc_df.loc[:, column].argsort().to_numpy()
    np.save(f"data/sorted/{column}_smiles_reorder.npy",
            smiles_reorder)

    sorted_vals = all_zinc_df.loc[:, column].iloc[smiles_reorder]
    np.save(f"data/sorted/{column}_values.npy",
            sorted_vals)

