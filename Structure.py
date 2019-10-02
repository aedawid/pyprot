
class Structure(object):

    def __init__(self, structure = None):
        


class Chain(object):

    def __init__(self, seq, bond_length = 0.38):
        self.beads         = [AA(aa) for aa in seq]
        self.name          = [resid.name       for resid in self.beads]
        self.id            = [resid.id         for resid in self.beads]
        self.weight        = [resid.weight     for resid in self.beads]
        self.frequency     = [resid.frequency  for resid in self.beads]
        self.charge        = [resid.charge     for resid in self.beads]
        self.polarity      = [resid.polarity   for resid in self.beads]
        self.aromatic      = [resid.aromatic   for resid in self.beads]
        self.hp_KD         = [resid.hp_KD      for resid in self.beads]

        self.seq_length    = len(seq)
        self.nbonds        = len(seq)-1
        self.bond_length   = bond_length


class Residue(object):

    def __init__(self, aa):
        if len(aa) == 3:
            aa = aa_names[resid]
        self.name = aa
        self.id = aa_features[aa][0]
        self.weight = aa_features[aa][1]
        self.frequency = aa_features[aa][2]
        self.charge = aa_features[aa][3]
        self.polarity = aa_features[aa][4]
        self.aromatic = aa_features[aa][5]
        self.hp_KD = aa_features[aa][6]

        self.x = 0.0
        self.y = 0.0
        self.z = 0.0


class Atom(object):

    def __init__(self, a):
        
