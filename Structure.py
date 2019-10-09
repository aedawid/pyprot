import re
import numpy as np

from math import sqrt
from copy import deepcopy
from itertools import combinations
from string import ascii_uppercase
from collections import OrderedDict

from PDB_menager import AA_code

class Structure(object):

    def __init__(self, PDB = None, m_id = 0):
        """Creates a Structure object for chosen model from PDB object."""

        self.structure = []
        if PDB:
            if isinstance(PDB, list):
                self.structure.append(Atom(atom, num)) for num, atom in enumerate(PDB[m_id])
            elif isinstance(PDB, str):
                ch = PDB.split('/')
                for c in ch:
                    names = PDB.split(':')
                    if len(names) > 1:
                        self.structure.append(Atom(atom, num, names[0])) for num, atom in enumerate(names[1])
                    else:
                        self.structure.append(Atom(atom, num)) for num, atom in enumerate(names)

        self.chains = self._get_chains()
        self.residues = self._get_resids()

    def _get_chains(self):
        

    def _get_resids(self):
        

class Chain(object):

    def __init__(self, chain):
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


class UnitedAtom(object):

    def __init__(self, atoms, frame=0, shift=1):
        """Creates UnitedAtom from given list of atoms."""
        for idx, atom in enumerate(atoms):
            if idx % shift == 0:
                
        
        
        
        
        self.atom_id = self.resi_id = num
        self.atom_name = "CA"
        self.resi_name = AA_code(residue)
        self.chain_id = chain
        self.coordin = Vector3d(0.0, 0.0, 0.0)
        self.occupancy = 0.0
        self.bfactor = 0.0
        self.atom_type = 'C'
        
        
        self.weight = aa_features[aa][1]
        self.frequency = aa_features[aa][2]
        self.charge = aa_features[aa][3]
        self.polarity = aa_features[aa][4]
        self.aromatic = aa_features[aa][5]
        self.hp_KD = aa_features[aa][6]



class Atom(object):

    def __init__(self, atom=None, num=0, chain='A'):

        self.idx = self.atom_id = self.resi_id = num
        self.atom_name = "CA"
        self.resi_name = "NAN"
        self.chain_id = chain
        self.coordin = Vector3d(0.0, 0.0, 0.0)
        self.occupancy = 0.0
        self.bfactor = 0.0
        self.atom_type = 'C'
        if isinstance(atom, list):
            """Creates Atom from records in PDB object."""
            self.atom_id = int(atom[1])
            self.atom_name = atom[2].strip()
            self.resi_name = atom[3].strip()
            self.chain_id = atom[4]
            self.resi_id = int(atom[5])
            self.coordin = Vector3d(atom[6], atom[7], atom[8])
            self.occupancy = float(atom[9])
            self.bfactor = float(atom[10])
            self.atom_type = atom[11].replace('\n','')
        elif isinstance(atom, str):
            if len(atom) > 1:
                """Creates Atom from string in pdb format."""
                self.atom_id = int(atom[6:11])
                self.atom_name = atom[11:16].strip()
                self.resi_name = atom[17:21].strip()
                self.chain_id = atom[21]
                self.resi_id = int(atom[22:26])
                self.coordin = Vector3d(atom[30:38], atom[38:46], atom[46:54])
                self.occupancy = float(atom[54:60])
                self.bfactor = float(atom[60:66])
                self.atom_type = atom[66:].replace('\n','')
            else:
                """Creates Atom from amino acid sequence"""
                self.resi_name = AA_code(residue)

#        self.aa_code = AA_code(self.resi_name)
#        self.weight = AA_ATTRIBUTES[self.aa_code][1]
#        self.frequency = AA_ATTRIBUTES[self.aa_code][2]
#        self.charge = AA_ATTRIBUTES[self.aa_code][3]
#        self.polarity = AA_ATTRIBUTES[self.aa_code][4]
#        self.aromatic = AA_ATTRIBUTES[self.aa_code][5]
#        self.hp_KD = AA_ATTRIBUTES[self.aa_code][6]

    def __str__(self):
        line = "ATOM  "
        name = " %-3s" % self.atom_name
        if len(self.name) == 4:
            name = self.atom_name
        line += "%5d %4s %-4s%1s%4d    %24s%6.2f%6.2f %s" % (
                self.serial, name, self.resi_name, self.chain_id, self.resi_id,
                self.coordin, self.occupancy, self.bfactor, self.atom_type
        )
        return line

    def is_chain(self, atom):
        return self.chain_id == atom.chain_id

    def is_resi(self, atom):
        return self.is_chain(atom) and self.resi_id == atom.resi_id

    def dist2(self, atom):
        return (self.coordin - atom.coordin).mod2()

    def distance(self, atom):
        return sqrt(self.dist2(atom))

    def min_dist(self, atoms):
        return min(self.distance(atom) for atom in atoms)

    def max_dist(self, atoms):
        return max(self.distance(atom) for atom in atoms)

    def add_coordin(self, atom):
        if isinstance(atom, Vector3d):
            self.coordin += atom
        else:
            self.coordin += atom.coordin
