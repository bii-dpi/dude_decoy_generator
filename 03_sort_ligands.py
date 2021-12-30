"""Write pre-sorted ZINC15 compound lists by each of the six properties."""

import numpy as np
import pandas as pd
from os import makedirs
from shutil import rmtree
from progressbar import progressbar


rmtree("data/sorted", ignore_errors=True)
makedirs("data/sorted")


print(f"Sorting SMILES by the relevant properties to make property-matching of"
      f" decoys faster...")

COLUMNS = ("smiles", "mol_wt", "logp",
           "num_rotatable", "num_hba", "num_hbd", "net_charge")

all_zinc_df = pd.read_csv("data/all_zinc_w_prop.csv",
                          names=COLUMNS,
                          dtype={"smiles": str, "mol_wt": np.float32,
                                 "logp": np.float32, "num_rotatable": np.int16,
                                 "num_hba": np.int16, "num_hbd": np.int16,
                                 "net_charge": np.int16})

# Save the reference ordering of the SMILES.
np.save("data/all_zinc_smiles_ref.npy", all_zinc_df[["smiles"]].to_numpy())

# For each of the six properties, save the relative ordering of the SMILES and
# property value list after sorting these values in ascending order.
for column in progressbar(COLUMNS[1:]):
    smiles_reorder = all_zinc_df.loc[:, column].argsort().to_numpy().flatten()
    np.save(f"data/sorted/{column}_smiles_reorder.npy",
            smiles_reorder)

    sorted_vals = all_zinc_df.loc[:, column].iloc[smiles_reorder].to_numpy().flatten()
    np.save(f"data/sorted/{column}_values.npy",
            sorted_vals)

print("Done.")

