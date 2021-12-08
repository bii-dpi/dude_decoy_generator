from rdkit import Chem
from rdkit.Chem import rdPartialCharges

a = rdPartialCharges.ComputeGasteigerCharges(Chem.MolFromSmiles("Cc1ccccc1",
sanitize=True))
print(type(a))
