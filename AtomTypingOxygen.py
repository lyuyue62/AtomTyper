#!/usr/bin/env python
import sys
sys.dont_write_bytecode = True

from Edge import Edge
from Atom import Atom
from Bond import Bond

class AtomTypingOxygen(object):
    def __init__(self):
        """__init__"""

    def setAtomTypeForOxygens(self, mol):
        num_atoms = len(mol.atoms)
        atom = Atom()
        atom_linked = Atom()
        index = -1
        i = 0
        while i < num_atoms:
            atom = mol.atoms[i]
            if atom.element.lower() == "O".lower():
                atom_linked = mol.atoms[atom.linkage[0]]
                if not atom.isRingAtom:
                    if atom_linked.isCarbonyl:
                        if atom_linked.numRingAtoms[0] == 6 and atom_linked.isAromatic:
                            #  OG2D4: 6-mem aromatic carbonyl oxygen (nucleic bases)
                            atom.atomType = "OG2D4"
                        elif atom_linked.numRingAtoms[0] == 5:
                            if not atom_linked.isKetone and (atom_linked.isAmide or atom_linked.isEster):
                                if not atom.isDeprotonatedOxygen:
                                    #  OG2D1: carbonyl O: amides, esters, [neutral] carboxylic acids, aldehydes, uera
                                    atom.atomType = "OG2D1"
                        elif atom.num_linkages == 2 and atom.numHydrogenAtoms == 1:
                            #  OG311: hydroxyl oxygen
                            atom.atomType = "OG311"
                        else:
                            if not atom_linked.isKetone and (atom_linked.isAmide or atom_linked.isEster):
                                if not atom.isDeprotonatedOxygen:
                                    #  OG2D1: carbonyl O: amides, esters, [neutral] carboxylic acids, aldehydes, uera
                                    atom.atomType = "OG2D1"
                                else:
                                    #  OG2D2: carbonyl O: negative groups: carboxylates, carbonate
                                    atom.atomType = "OG2D2"
                            else:
                                #  OG2D3: carbonyl O: ketones
                                atom.atomType = "OG2D3"
                    elif self.isTypeOG2D5(mol, i):
                        #  OG2D5: CO2 oxygen
                        atom.atomType = "OG2D5"
                    elif atom_linked.atomType != None and atom_linked.atomType.lower() == "NG2O1".lower():
                        #  OG2N1: NITB, nitrobenzene
                        atom.atomType = "OG2N1"
                    elif atom.isPOxide or atom.isSOxide:
                        if atom.num_linkages == 1:
                            #  OG2P1: =O in phosphate or sulfate
                            atom.atomType = "OG2P1"
                        elif atom.num_linkages == 2 and atom.numHydrogenAtoms == 0:
                            #  OG303: phosphate/sulfate ester oxygen
                            atom.atomType = "OG303"
                        elif atom.num_linkages == 2 and atom.numHydrogenAtoms == 1:
                            #  OG311: hydroxyl oxygen
                            atom.atomType = "OG311"
                    elif atom.num_linkages == 2 and atom.numHydrogenAtoms == 1:
                        #  OG311: hydroxyl oxygen
                        atom.atomType = "OG311"
                    elif atom.num_linkages == 2 and atom.numCarbonAtoms == 2:
                        #  OG301: ether -O- !SHOULD WE HAVE A SEPARATE ENOL ETHER??? IF YES, SHOULD WE MERGE IT WITH OG3R60???
                        atom.atomType = "OG301"
                    if atom.isEster and atom.num_linkages == 2 and atom.numCarbonAtoms == 2:
                        #  OG302: ester -O-
                        atom.atomType = "OG302"
                    elif self.isTypeOG304(mol, i):
                        #  OG304: linkage oxygen in pyrophosphate/pyrosulphate
                        atom.atomType = "OG304"
                    elif not atom.isEster and atom.isDeprotonatedOxygen and (atom.numCarbonAtoms == 1 or atom.numNitrogenAtoms == 1):
                        #  OG312: ionized alcohol oxygen
                        atom.atomType = "OG312"
                else:
                    if atom.isFuran:
                        #  OG2R50: FURA, furan
                        atom.atomType = "OG2R50"
                    elif atom.isFuranose:
                        #  OG3C51: 5-mem furanose ring oxygen (ether)
                        atom.atomType = "OG3C51"
                    elif atom.numRingAtoms[0] == 6:
                        #  OG3R60: O in 6-mem cyclic enol ether (PY01, PY02) or ester
                        atom.atomType = "OG3R60"
                        if not self.linkedCarbonHasDoundBond(mol, i):
                            #  OG3C61: DIOX, dioxane, ether in 6-membered ring !SHOULD WE MERGE THIS WITH OG3R60???
                            atom.atomType = "OG3C61"
                    elif atom.numRingAtoms[0] == 3 and atom.num_linkages == 2 and atom.numCarbonAtoms == 2:
                        #  OG3C31: epoxide oxygen, 1EOX, 1BOX, sc
                        atom.atomType = "OG3C31"
            i += 1

    def getIndexAtomType(self, mol, atm_index, atomType):
        atom = mol.atoms[atm_index]
        i = 0
        while i < 5:
            if atom.linkage[i] == -1:
                return -1
            else:
                if mol.atoms[atom.linkage[i]].atomType != None and mol.atoms[atom.linkage[i]].atomType.lower() == atomType.lower():
                    return atom.linkage[i]
            i += 1
        return -1

    def getIndexElement(self, mol, atm_index, element):
        atom = mol.atoms[atm_index]
        i = 0
        while i < 5:
            if atom.linkage[i] == -1:
                return -1
            else:
                if mol.atoms[atom.linkage[i]].element.lower() == element.lower():
                    return atom.linkage[i]
            i += 1
        return -1

    def isTypeOG2D5(self, mol, atm_index):
        atom = mol.atoms[atm_index]
        atom_linked = Atom()
        if atom.num_linkages == 1:
            atom_linked = mol.atoms[atom.linkage[0]]
            if atom_linked.element.lower() == "C".lower() and atom_linked.num_linkages == 2 and atom_linked.numOxygenAtoms == 2:
                return True
            elif atom_linked.element.lower() == "C".lower() and atom_linked.num_linkages == 2 and atom_linked.numOxygenAtoms == 1 and atom_linked.numNitrogenAtoms == 1:
                return True
        return False

    def isTypeOG304(self, mol, atm_index):
        atom = mol.atoms[atm_index]
        atom_linked_1 = Atom()
        atom_linked_2 = Atom()
        if atom.num_linkages == 2:
            atom_linked_1 = mol.atoms[atom.linkage[0]]
            atom_linked_2 = mol.atoms[atom.linkage[1]]
            if (atom_linked_1.element.lower() == "P".lower() or atom_linked_1.element.lower() == "S".lower()) and atom_linked_1.numOxygenAtoms == 4:
                if (atom_linked_2.element.lower() == "P".lower() or atom_linked_2.element.lower() == "S".lower()) and atom_linked_2.numOxygenAtoms == 4:
                    return True
        return False

    def linkedCarbonHasDoundBond(self, mol, atm_index):
        atom = mol.atoms[atm_index]
        if atom.num_linkages == 2 and atom.numCarbonAtoms == 2:
            if self.hasDoubleBond(mol, atom.linkage[0]) or self.hasDoubleBond(mol, atom.linkage[1]):
                return True
        return False

    def hasDoubleBond(self, mol, atm_index):
        atom = mol.atoms[atm_index]
        i = 0
        while i < atom.num_linkages:
            if atom.bondOrder[i] == 2 or atom.bondType[i].lower() == "ar".lower():
                return True
            i += 1
        return False

