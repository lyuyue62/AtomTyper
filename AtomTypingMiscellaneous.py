#!/usr/bin/env python
import sys
sys.dont_write_bytecode = True

from Edge import Edge
from Atom import Atom
from Bond import Bond

class AtomTypingMiscellaneous(object):
    def __init__(self):
        """__init__"""
    def setAtomTypeForMiscellaneousAtoms(self, mol):
        num_atoms = len(mol.atoms)
        atom = Atom()
        atom_linked = Atom()
        index = -1
        i = 0
        while i < num_atoms:
            atom = mol.atoms[i]
            if atom.element.lower() == "P".lower():
                if atom.numSulfurAtoms == 0:
                    charge = self.getCharge( mol, i )
                    if charge == 0:
                        #  PG0: neutral phosphate
                        atom.atomType = "PG0"
                    elif charge == -1:
                        #  PG1: phosphate -1
                        atom.atomType = "PG1"
                    elif charge == -2:
                        #  PG2: phosphate -2
                        atom.atomType = "PG2"
                else:
                    j = 0
                    while j < atom.num_linkages:
                        atom_linked = mol.atoms[atom.linkage[j]]
                        if atom_linked.element.lower() == "O".lower() and atom_linked.num_linkages == 1:
                            #  OG2S1: mono-thio S-P bond modulated oxygen; lsk
                            atom_linked.atomType = "OG2S1"
                        j += 1
            elif atom.element.lower() == "LP".lower():
                #  LPH: Lone pair for halogens
                atom.atomType = "LPH"
            elif atom.element.lower() == "AL".lower():
                if self.countSpecificElement(mol, i, "F") == 4:
                    #  ALG1: Aluminum, for ALF4, AlF4-
                    atom.atomType = "ALG1"
            elif atom.element.lower() == "S".lower():
                if atom.num_linkages == 1:
                    atom_linked = mol.atoms[atom.linkage[0]]
                    if atom_linked.element.lower() == "P".lower():
                        if num_S == 1:
                            #  SG2P1: mono-thio S-P bond; lsk
                            atom.atomType = "SG2P1"
                        elif num_S == 2:
                            #  SG2P2: di-thio S-P bond; lsk
                            atom.atomType = "SG2P2"
            i += 1

    def getIndex(self, mol, atm_index, element):
        atom = mol.atoms[atm_index]
        i = 0
        while i < atom.num_linkages:
            if mol.atoms[atom.linkage[i]].element.lower() == element.lower():
                return atom.linkage[i]
            i += 1
        return -1

    def countSpecificElement(self, mol, atm_index, element):
        atom = mol.atoms[atm_index]
        num_elements = 0
        i = 0
        while i < atom.num_linkages:
            if mol.atoms[atom.linkage[i]].element.lower() == element.lower():
                num_elements += 1
            i += 1
        return num_elements

    def getCharge(self, mol, atm_index):
        atom = mol.atoms[atm_index]
        atom_linked = Atom()
        num_end_Oxygens = 0
        num_end_Oxygens_w_doubleBond = 0
        i = 0
        while i < atom.num_linkages:
            atom_linked = mol.atoms[atom.linkage[i]]
            if atom_linked.num_linkages == 1:
                num_end_Oxygens += 1
                if atom_linked.bondOrder[0] == 2:
                    num_end_Oxygens_w_doubleBond += 1
            i += 1
        return num_end_Oxygens_w_doubleBond - num_end_Oxygens

