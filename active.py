"""Lorem ipsum."""

class Active():
    def __init__(self, pair, assignable_candidates):
        self.pdb_id, self.smiles = pair
        self.assignable_candidates = {assignable_smiles: None
                                      for assignable_smiles in
                                      assignable_candidates}
        self.assigned_decoys = []


    def has_50_decoys(self):
        return len(self) == 50


    def can_assign(self, cand_smiles):
        return cand_smiles in self.assignable_candidates


    def add_decoy(self, cand_decoy):
        self.assigned_decoys.append(cand_decoy)
        del self.assignable_candidates[cand_decoy]


    def get_rows(self):
        rows = []
        for assigned in self.assigned_decoys:
            rows.append(f"{assigned} {self.smiles} {self.pdb_id}")
        return rows


    def __len__(self):
        return len(self.assigned_decoys)


    def __eq__(self, other):
        return len(self) == len(other) and self.smiles == other.smiles


    def __ne__(self, other):
        return len(self) != len(other)


    def __lt__(self, other):
        if len(self) == len(other):
            return (len(self.assignable_candidates) <
                    len(other.assignable_candidates))
        return len(self) < len(other)


    def __le__(self, other):
        if len(self) == len(other):
            return (len(self.assignable_candidates) <=
                    len(other.assignable_candidates))
        return len(self) <= len(other)


    def __gt__(self, other):
        if len(self) == len(other):
            return (len(self.assignable_candidates) >
                    len(other.assignable_candidates))
        return len(self) > len(other)


    def __ge__(self, other):
        if len(self) == len(other):
            return (len(self.assignable_candidates) >=
                    len(other.assignable_candidates))
        return len(self) >= len(other)

