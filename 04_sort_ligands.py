import numpy as np

import pandas as pd
from os import remove
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor

all_zinc_df = pd.read_csv("data/all_zinc_w_prop.csv",
                          names=("smiles", "mol_wt", "logp",
                                       "num_rotatable", "num_hba", "num_hbd"),
                          dtype={"smiles": str, "mol_wt": np.float32,
                                 "logp": np.float32, "num_rotatable": int,
                                 "num_hba": int, "num_hbd": int})

try:
    os.remove(r"data/all_zinc_w_prop.db")
except:
    pass


print("Done")

'''
cursor = conn.cursor()
cursor.execute("""CREATE TABLE all_prop
               (smiles TEXT,
                mol_wt REAL,
                logp INT,
                num_rotatable INT,
                num_hba INT,
                num_hbd INT);""")

with open("data/all_zinc_w_prop.smi","r") as f:
    dr = csv.DictReader(f, fieldnames=("smiles", "mol_wt", "logp",
                                       "num_rotatable", "num_hba", "num_hbd"))
    to_db = [(i["smiles"], i["mol_wt"], i["logp"],
              i["num_rotatable"], i["num_hba"], i["num_hbd"]) for i in dr]

cursor.executemany("""INSERT INTO all_prop
                      (smiles, mol_wt, logp, num_rotatable, num_hba, num_hbd)
                      VALUES (?, ?, ?, ?, ?, ?);""",
                   to_db)
conn.commit()
conn.close()
'''

'''
ExactMolWt(mol), MolLogP(mol),
 CalcNumRotatableBonds(mol),
 CalcNumHBA(mol), CalcNumHBD(mol)

class Flexlist(list):
    def __getitem__(self, keys):
        if isinstance(keys, (int, slice)): return list.__getitem__(self, keys)
        return Flexlist([[self[i][k] for k in keys]
                         for i in range(len(self))])


class Sort(object):
    def sort(self, l):
        l.sort()
        return l

    def isSorted(self, a):
        for i in range(1, len(a)):
            if a[i-1] > a[i]:
                return False
        return True

    def time(self, l):
        def real():
            self.sort(l)
        return timeit.timeit(real, number =1)


class QuickSort(Sort):
    def sorthelper(self, a , s, e):
        if (e-s)==0:
            return
        pivot = a[e-1]
        p1 = s
        p2 = e - 1
        while (p1 != p2):
            if (a[p1] > pivot):
                a[p2] = a[p1]
                a[p1] = a[p2-1]
                a[p2-1] = pivot
                p2 = p2 -1
            else:
                p1+=1
        self.sorthelper(a, s, p2)
        self.sorthelper(a, p2+1, e)

    def sort(self, a):
        self.sorthelper(a , 0, len(a))


if __name__ == "__main__":
    def write_sorted(pair):
        i, zinc_list = pair
        print("Sorting")
        sorted_list = [" ".join([str(x) for x in zinc_list[i]])
                       for i in np.argsort(np.array(zinc_list[[-1]],
                                                    dtype=np.float16))]
        print("Writing")
        with open(f"data/sorted_zinc_{i}", "w") as f:
            f.write("\n".join(sorted_list))


    print("Sorting and saving all ZINC compounds according to properties.")

    print("Loading ZINC compound list...")
    with open("data/all_zinc_w_prop.smi", "r") as f:
        all_zinc_w_prop = [line.strip("\n").split() for line in f.readlines()]
    all_zinc_w_prop = [[line[0]] + list(map(lambda x: float(x), line[1:]))
                       for line in all_zinc_w_prop]
    all_zinc_w_prop = Flexlist(all_zinc_w_prop)

    zinc_copies = [(i, all_zinc_w_prop[[0, i]]) for i in range(1, 6)]

    print("Saving sorted ZINC compound lists...")
    """
    with ProcessPoolExecutor() as executor:
        executor.map(write_sorted, zinc_copies)
    """
    for pair in tqdm(zinc_copies):
        write_sorted(pair)

    print("Done.")

'''
