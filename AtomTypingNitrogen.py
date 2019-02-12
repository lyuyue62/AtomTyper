#!/usr/bin/env python
import sys
sys.dont_write_bytecode = True

from Edge import Edge
from Atom import Atom
from Bond import Bond

class AtomTypingNitrogen(object):
    def __init__(self):
        """__init__"""

    def setAtomTypingNitrogen(self, mol):
        num_atoms = len(mol.atoms)
        atom = Atom()
        iminCarbonIndex = -1
        i = 0
        while i < num_atoms:
            atom = mol.atoms[i]
            if atom.element.lower() == "N".lower():
                #  NG1T1: N for cyano group
                if self.getIndexAtomType(mol, i, "CG1N1") != -1:
                    atom.atomType = "NG1T1"
                #  NG2D1: N for neutral imine/Schiff's base (C=N-R, acyclic amidine, gunaidine)
                #  NG2P1: N for protonated imine/Schiff's base (C=N(+)H-R, acyclic amidinium, guanidinium)
                iminCarbonIndex = self.getIndexOfImineCarbon(mol, i)
                if iminCarbonIndex != -1:
                    if not mol.atoms[iminCarbonIndex].isProtonatedImineGroup:
                        atom.atomType = "NG2D1"
                    else:
                        atom.atomType = "NG2P1"
                if atom.isAmide:
                    if atom.num_linkages == 3:
                        if atom.numHydrogenAtoms == 0:
                            #  NG2S0: N,N-disubstituted amide, proline N (CO=NRR')
                            atom.atomType = "NG2S0"
                        elif atom.numHydrogenAtoms == 1:
                            #  NG2S1: peptide nitrogen (CO=NHR)
                            atom.atomType = "NG2S1"
                        elif atom.numHydrogenAtoms == 2:
                            #  NG2S2: terminal amide nitrogen (CO=NH2)
                            atom.atomType = "NG2S2"
                if atom.isRingAtom:
                    if not atom.isBridgingAtom:
                        if atom.isAmide and atom.numRingAtoms[0] == 4:
                            #  NG2R43: amide in 4-memebered ring (planar), AZDO, lsk
                            atom.atomType = "NG2R43"
                        elif atom.numRingAtoms[0] == 5:
                            if not atom.isProtonatedNitrogen and atom.num_linkages == 2 and atom.numHydrogenAtoms == 0:
                                #  NG2R50: double bound neutral 5-mem planar ring, purine N7
                                atom.atomType = "NG2R50"
                            elif not atom.isProtonatedNitrogen and atom.num_linkages == 3:
                                if atom.isAromatic:
                                    #  NG2R51: single bound neutral 5-mem planar (all atom types sp2) ring, his, trp pyrrole (fused)
                                    atom.atomType = "NG2R51"
                                if isTypeNG2R53(mol, i):
                                    #  NG2R53: amide in 5-memebered NON-SP2 ring (slightly pyramidized), 2PDO, kevo
                                    atom.atomType = "NG2R53"
                                if isBipyrroleNitrogen(mol, i):
                                    #  NG2R57: 5-mem ring, bipyrroles
                                    atom.atomType = "NG2R57"
                                if not atom.isAromatic and (atom.numHydrogenAtoms == 0 or atom.numHydrogenAtoms == 1) and not atom.isAmide:
                                    #  NG3C51: secondary sp3 amine in 5-membered ring
                                    atom.atomType = "NG3C51"
                            elif iminCarbonIndex != -1 and mol.atoms[iminCarbonIndex].isProtonatedImineGroup:
                                #  NG2R52: protonated schiff base, amidinium, guanidinium in 5-membered ring, HIS, 2HPP, kevo
                                atom.atomType = "NG2R52"
                            elif atom.isAromatic and atom.num_linkages == 3 and not hasDoubleBond(mol, i):
                                #  NG2R51: single bound neutral 5-mem planar (all atom types sp2) ring, his, trp pyrrole (fused)
                                atom.atomType = "NG2R51"
                        elif atom.numRingAtoms[0] == 6:
                            if atom.isAromatic and atom.numHydrogenAtoms == 0:
                                #  NG2R60: double bound neutral 6-mem planar ring, pyr1, pyzn
                                atom.atomType = "NG2R60"
                            if atom.isAromatic and (hasCarbonyl(mol, i) or atom.isProtonatedNitrogen):
                                #  NG2R61: single bound neutral 6-mem planar ring imino nitrogen; glycosyl linkage
                                atom.atomType = "NG2R61"
                            if isTypeNG2R62(mol, i):
                                #  NG2R62: double bound 6-mem planar ring with heteroatoms in o or m, pyrd, pyrm
                                atom.atomType = "NG2R62"
                            if isTypeNG2R67(mol, i):
                                #  NG2R67: 6-mem planar ring substituted with 6-mem planar ring (N-phenyl pyridinones etc.)
                                atom.atomType = "NG2R67"
                    else:
                        if isTypeNG2RC0(mol, i):
                            #  NG2RC0: 6/5-mem ring bridging N, indolizine, INDZ, kevo
                            atom.atomType = "NG2RC0"
                if not atom.isAmide and not atom.isProtonatedNitrogen and iminCarbonIndex == -1:
                    if isTypeNG301(mol, i):
                        #  NG301: neutral trimethylamine nitrogen
                        atom.atomType = "NG301"
                    elif atom.numRingAtoms[0] != 5 and isTypeNG311(mol, i):
                        #  NG311: neutral dimethylamine nitrogen
                        atom.atomType = "NG311"
                    elif isTypeNG321(mol, i):
                        #  NG321: neutral methylamine nitrogen
                        atom.atomType = "NG321"
                    elif atom.num_linkages == 3 and atom.numHydrogenAtoms == 3:
                        #  NG331: neutral ammonia nitrogen
                        atom.atomType = "NG331"
                    elif isTypeNG3N1(mol, i):
                        #  NG3N1: N in hydrazine, HDZN
                        atom.atomType = "NG3N1"
                #  NG2S3: external amine ring nitrogen (planar/aniline), phosphoramidate (PO3-NR2), sulfamate (SO3-NR2)
                if self.isTypeNG2S3(mol, i):
                    atom.atomType = "NG2S3"
                #  NG2O1: NITB, nitrobenzene
                if atom.num_linkages == 3 and atom.numOxygenAtoms == 2 and atom.numCarbonAtoms == 1:
                    atom.atomType = "NG2O1"
                if not atom.isRingAtom and atom.isProtonatedNitrogen and atom.numHydrogenAtoms == 0:
                    #  NG3P0: quarternary N+, choline
                    atom.atomType = "NG3P0"
                elif atom.isRingAtom and not atom.isAromatic and atom.isProtonatedNitrogen and atom.numHydrogenAtoms == 1:
                    #  NG3P1: tertiary NH+ (PIP)
                    atom.atomType = "NG3P1"
                elif atom.isRingAtom and not atom.isAromatic and atom.isProtonatedNitrogen and atom.numHydrogenAtoms == 2:
                    #  NG3P2: secondary NH2+ (proline)
                    atom.atomType = "NG3P2"
                elif not atom.isRingAtom and atom.isProtonatedNitrogen and atom.numHydrogenAtoms == 3:
                    #  NG3P3: primary NH3+, phosphatidylethanolamine
                    atom.atomType = "NG3P3"
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

    def countSpecificElement(self, mol, atm_index, element):
        atom = mol.atoms[atm_index]
        num_elements = 0
        i = 0
        while i < 5:
            if atom.linkage[i] != -1 and mol.atoms[atom.linkage[i]].element.lower() == element.lower():
                num_elements += 1
            i += 1
        return num_elements

    def getIndexDoubleBondedCarbon(self, mol, atm_index):
        atom = mol.atoms[atm_index]
        i = 0
        while i < 4:
            if atom.linkage[i] != -1 and atom.bondOrder[i] == 2 and mol.atoms[atom.linkage[i]].element.lower() == "C".lower():
                return atom.linkage[i]
            i += 1
        return -1

    def hasDifferentRingIndex(self, mol, atom_i, atom_j):
        index_i = int()
        index_j = int()
        i = 0
        while i < 3:
            index_i = mol.atoms[atom_i].ringIndex[i]
            while j < 3:
                index_j = mol.atoms[atom_j].ringIndex[j]
                if index_i != -1 and index_j != -1 and index_i != index_j:
                    return True
                j += 1
            i += 1
        return False

    def hasSameRingIndex(self, mol, atom_i, atom_j):
        index_i = int()
        index_j = int()
        i = 0
        while i < 3:
            index_i = mol.atoms[atom_i].ringIndex[i]
            while j < 3:
                index_j = mol.atoms[atom_j].ringIndex[j]
                if index_i != -1 and index_j != -1 and index_i == index_j:
                    return True
                j += 1
            i += 1
        return False

    def hasDoubleBond(self, mol, atm_index):
        atom = mol.atoms[atm_index]
        i = 0
        while i < 4:
            if atom.linkage[i] != -1 and atom.bondOrder[i] == 2:
                return True
            i += 1
        return False

    def getIndexOfImineCarbon(self, mol, atm_index):
        atom = mol.atoms[atm_index]
        atom_linked = Atom()
        if atom.isRingAtom:
            i = 0
            while i < atom.num_linkages:
                atom_linked = mol.atoms[atom.linkage[i]]
                if atom_linked.isImineCarbon and atom_linked.isRingAtom:
                    return atom.linkage[i]
                i += 1
        else:
            i = 0
            while i < atom.num_linkages:
                atom_linked = mol.atoms[atom.linkage[i]]
                if atom_linked.isImineCarbon and not atom_linked.isRingAtom:
                    return atom.linkage[i]
                i += 1
        return -1

    def isBipyrroleNitrogen(self, mol, atm_index):
        atom = mol.atoms[atm_index]
        atom_linked = Atom()
        if atom.isPyrrole:
            while i < atom.num_linkages:
                atom_linked = mol.atoms[atom.linkage[i]]
                if atom_linked.isPyrrole and not mol.atoms[atom.linkage[i]].isBridgingAtom and atom.ringIndex[0] != atom_linked.ringIndex[0]:
                    return True
                i += 1
        return False

    def hasCarbonyl(self, mol, atm_index):
        atom = mol.atoms[atm_index]
        i = 0
        while i < 3:
            if atom.linkage[i] != -1 and mol.atoms[atom.linkage[i]].isCarbonyl:
                return True
            i += 1
        return False

    def isTypeNG2S3(self, mol, atm_index):
        atom = mol.atoms[atm_index]
        atom_linked = Atom()
        if not atom.isRingAtom:
            if atom.num_linkages == 3:
                i = 0
                while i < 3:
                    atom_linked = mol.atoms[atom.linkage[i]]
                    if atom_linked.element.lower() == "C".lower() and atom_linked.isRingAtom and atom_linked.isAromatic and atom.numHydrogenAtoms == 2:
                        return True
                    elif (atom_linked.element.lower() == "P".lower() or atom_linked.element.lower() == "S".lower()) and self.countSpecificElement(mol, atom.linkage[i], "O") == 3:
                        return True
                    i += 1
        return False

    def isTypeNG2R62(self, mol, atm_index):
        atom = mol.atoms[atm_index]
        atom_linked = Atom()
        if atom.isAromatic and atom.num_linkages == 2:
            while i < 2:
                atom_linked = mol.atoms[atom.linkage[i]]
                if atom_linked.isAromatic:
                    if not atom_linked.element.lower() == "C".lower():
                        return True
                    else:
                        while j < atom_linked.num_linkages:
                            if atom_linked.linkage[j] != atm_index and mol.atoms[atom_linked.linkage[j]].isAromatic and not mol.atoms[atom_linked.linkage[j]].element.lower() == "C".lower():
                                return True
                            j += 1
                i += 1
        return False

    def isTypeNG2R67(self, mol, atm_index):
        atom = mol.atoms[atm_index]
        index_i = -1
        index_j = -1
        i = 0
        while i < atom.num_linkages:
            if atom.numRingAtoms[0] == 6 and mol.atoms[atom.linkage[i]].numRingAtoms[0] == 6 and not mol.atoms[atom.linkage[i]].isBridgingAtom and atom.ringIndex[0] != mol.atoms[atom.linkage[i]].ringIndex[0]:
                return True
            i += 1
        return False

    def isTypeNG2RC0(self, mol, atm_index):
        atom = mol.atoms[atm_index]
        isFiveMemberedRingAtom = False
        isSixMemberedRingAtom = False
        if atom.isAromatic:
            while i < 3:
                if atom.numRingAtoms[i] == 5:
                    isFiveMemberedRingAtom = True
                if atom.numRingAtoms[i] == 6:
                    isSixMemberedRingAtom = True
                i += 1
            if isFiveMemberedRingAtom and isSixMemberedRingAtom:
                return True
        return False

    def isTypeNG301(self, mol, atm_index):
        atom = mol.atoms[atm_index]
        if not atom.isAromatic and (atom.numCarbonAtoms == 3 and atom.numHydrogenAtoms == 0):
            return True
        elif not atom.isAromatic and (atom.numCarbonAtoms == 2 and atom.numSulfurAtoms == 1 and atom.numHydrogenAtoms == 0):
            return True
        return False

    def isTypeNG311(self, mol, atm_index):
        atom = mol.atoms[atm_index]
        if atom.numCarbonAtoms == 2 and atom.numHydrogenAtoms == 1:
            return True
        elif atom.numCarbonAtoms == 1 and (atom.numSulfurAtoms == 1 or atom.numOxygenAtoms == 1) and atom.numHydrogenAtoms == 1:
            return True
        return False

    def isTypeNG321(self, mol, atm_index):
        atom = mol.atoms[atm_index]
        if (atom.numCarbonAtoms == 1 or atom.numSulfurAtoms == 1) and atom.numHydrogenAtoms == 2:
            return True
        return False

    def isTypeNG3N1(self, mol, atm_index):
        atom = mol.atoms[atm_index]
        atom_linked = Atom()
        if atom.num_linkages == 3:
            if index != -1:
                atom_linked = mol.atoms[index]
                if atom_linked.numHydrogenAtoms == 2 or atom_linked.numHydrogenAtoms == 1:
                    return True
        return False

    def isTypeNG2R53(self, mol, atm_index):
        atom = mol.atoms[atm_index]
        atom_linked = Atom()
        if atom.isAmide:
            while i < atom.num_linkages:
                atom_linked = mol.atoms[atom.linkage[i]]
                if atom_linked.isCarbonyl and atom.ringIndex[0] == atom_linked.ringIndex[0]:
                    return True
                i += 1
        return False

