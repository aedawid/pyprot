import re
import numpy as np
import hashlib

from math import sqrt
from copy import deepcopy
from itertools import combinations
from string import ascii_uppercase
from collections import OrderedDict

from PDB_menager import AA_code, AA_ATTRIBUTES
from vector3d import Vector3d

class Structure(object):

    def __init__(self, PDB = None, m_id = 0):
        """Creates a Structure object for chosen model from PDB object."""

        self.structure = []
        if PDB:
            if isinstance(PDB, list):
                if isinstance(PDB[0], Atom):
                    self.structure = PDB
                else:
                    for num, atom in enumerate(PDB[m_id]):
                        self.structure.append(Atom(atom, num))
            elif isinstance(PDB, str):
                ch = PDB.split('/')
                for c in ch:
                    names = PDB.split(':')
                    if len(names) > 1:
                        for num, atom in enumerate(names[1]):
                            self.structure.append(Atom(atom, num, names[0]))
                    else:
                        for num, atom in enumerate(names):
                            self.structure.append(Atom(atom, num))

        self.chains   = self._get_chains()
        self.residues = self._get_resids()
        self.sequence = self._get_seq()
        self.seq_len  = len(self.sequence)
        self.united_codes = AA_UNITED
        self.united_attributes = self.calculate_united_attributes()


    def _get_chains(self):
        """Returns dictionary of chains in Structure object and its atom ranges."""
        chains = {}
        for num, atom in enumerate(self.structure):
            if atom.chain_id in chains.keys():
                chains[atom.chain_id][1] = num
            else:
                chains.update({atom.chain_id : [num, 0]})
        return chains

    def _get_resids(self):
        """Returns dictionary of residues in Structure object and its atom ranges."""
        resids = OrderedDict()
        key =''
        for num, atom in enumerate(self.structure):
            prefix = str(atom.resi_id)+":"+atom.chain_id
            if prefix in resids.keys():
                resids[prefix][1] = num+1
            else:
                resids[prefix] = [num, 0]
        return resids

    def _get_seq(self):
        """Returns amino acid sequence."""
        seq = ''
        for resid in self.residues:
            seq += self.structure[self.residues[resid][0]].aa_code
        return seq

    def select_atoms(self, atoms):
        """Returns selected group of atoms by its chain, names or range."""
        subset=[]
        if isinstance(atoms, str):
            for ch in self.chains:
                if ch == atoms:
                    subset = structure[self.chains[ch][0]:self.chains[ch][1]]
        if isinstance(atoms, tuple):
            subset = structure[atom[0]:atom[1]]
        else:
            for atom in self.structure:
                if atom.atom_name in atoms:
                    subset.append(atom)
        return subset

    def create_united_representation(self, frame=2, shift=1):
        """Creates united representation (structure object) of chosen frame and shift per residue."""
        n = num = 0
        subset = self.select_atoms(['CA'])
        representation=[]
        for atom in subset:
            is_ok = False
            if num + frame > len(subset):
                break
            if num % shift == 0:
                united_at = Atom()
                name  = atom.aa_code
                united_at.coordin = atom.coordin
                for i in range(1, frame):
                    if atom.is_chain(subset[num+i]):
                        is_ok = True

                        name += subset[num+i].aa_code
                        united_at.add_coordin(subset[num+i])
                    else:
                        is_ok = False
                if is_ok:
                    united_at.idx = united_at.atom_id = united_at.resi_id = n
                    n += 1
                    united_at.aa_code = name
                    united_at.chain_id = atom.chain_id
                    united_at.coordin.x /= float(frame)
                    united_at.coordin.y /= float(frame)
                    united_at.coordin.z /= float(frame)
                    if not name in AA_UNITED:
                        ss = hashlib.sha224(name.encode('utf-8')).hexdigest()[-3:]
                        if not ss in AA_UNITED.values():
                            AA_UNITED[name] = ss
                            united_at.resi_name = ss
                        else:
                            ss = hashlib.sha224(name.encode('utf-8')).hexdigest()
                            is_done = False
                            for i in range(0, len(ss)-3):
                                sss = ss[:3]
                                if not sss in AA_UNITED.values():
                                    AA_UNITED[name] = sss
                                    united_at.resi_name = sss
                                    is_done = True
                                    break
                            if is_done == False:
                                print("ERROR: hash value repeated!")
                representation.append(united_at)
            num += 1
        return Structure(representation)

    def calculate_united_attributes(self):
        """"""
        if not len(AA_UNITED):
            print("")#("AA_UNITED is empty.")
        else:
            for num, key in enumerate(AA_UNITED):
                weight = freq = charge = polar = aromatic = hp_KD = 0.0
                for aa in key:
                    weight += float(AA_ATTRIBUTES[aa][1])
                    freq += float(AA_ATTRIBUTES[aa][2])
                    charge += float(AA_ATTRIBUTES[aa][3])
                    polar += float(AA_ATTRIBUTES[aa][4])
                    aromatic += float(AA_ATTRIBUTES[aa][5])
                    hp_KD += float(AA_ATTRIBUTES[aa][6])
                UN_ATTRIBUTES[key] = (num, weight, freq, charge, polar, aromatic, hp_KD)
        return UN_ATTRIBUTES


    def place_attribute_in_bfcol(self, which=6):
        for atom in self.structure:
            atom.bfactor = UN_ATTRIBUTES[atom.aa_code][which]

    def place_seq_feature_in_bfcol(self, feature):
        vec = []
        if len(self.residues) != len(feature):
#            print("ERROR: Number of residues (", len(self.residues), ") and feature size (", len(feature), ") are not equal!")
            n = 0
            for num, res in enumerate(self.sequence):
                if res == feature[n][0]:
                    vec.append(feature[n][1])
                    n += 1
                elif n+3 < len(feature) and num+3 < len(self.sequence) and self.sequence[num+1] == feature[n+1][0] and self.sequence[num+2] == feature[n+2][0] and self.sequence[num+3] == feature[n+3][0]:
                    vec.append(0.000)
                    n += 1
                else:
                    vec.append(0.000)
        else:
            for ix in feature:
                vec.append(ix[1])
        for num, res in enumerate(self.residues):
            for atom in self.structure[self.residues[res][0]:self.residues[res][1]]:
                atom.bfactor = float(vec[num])

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

        self.aa_code = AA_code(self.resi_name)



    def __str__(self):
        line = "ATOM  "
        name = " %-3s" % self.atom_name
        if len(self.atom_name) == 4:
            name = self.atom_name
        line += "%5d %4s %-4s%1s%4d    %24s%6.2f%6.2f %s" % (
                self.atom_id, name, self.resi_name, self.chain_id, self.resi_id,
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

AA_UNITED = {}
UN_ATTRIBUTES = {}