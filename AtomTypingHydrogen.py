#!/usr/bin/env python
class AtomTypingHydrogen(object):
    def __init__(self):
        """__init__"""
    def setAtomTypeForHydrogens(self, mol):
        num_atoms = len(mol.atoms)
        atom = Atom()
        atom_linked = Atom()
        i = 0
        while i < num_atoms:
            atom = mol.atoms.elementAt(i)
            if atom.element.lower() == "H".lower():
                atom_linked = mol.atoms.elementAt(atom.linkage[0])
                if atom_linked.isAromatic:
                    #  HGR61: aromatic H
                    atom.atomType = "HGR61"
                if element.lower() == "C".lower():
                    if not atom_linked.isAromatic:
                        if atom_linked.numHydrogenAtoms == 1:
                            #  HGA1: alphatic proton, CH
                            atom.atomType = "HGA1"
                        elif atom_linked.numHydrogenAtoms == 2:
                            #  HGA2: alphatic proton, CH2
                            atom.atomType = "HGA2"
                        elif atom_linked.numHydrogenAtoms == 3:
                            #  HGA3: alphatic proton, CH3
                            atom.atomType = "HGA3"
                        if hasCarbonDoubleBond(mol, atom_linked):
                            if atom_linked.numHydrogenAtoms == 1:
                                #  HGA4: alkene proton; RHC=
                                atom.atomType = "HGA4"
                            elif atom_linked.numHydrogenAtoms == 2:
                                #  HGA5: alkene proton; H2C=CR
                                atom.atomType = "HGA5"
                        if num_fluorides > 0:
                            if num_fluorides == 1:
                                #  HGA6: aliphatic H on fluorinated C, monofluoro
                                atom.atomType = "HGA6"
                            elif num_fluorides == 2:
                                #  HGA7: aliphatic H on fluorinated C, difluoro
                                atom.atomType = "HGA7"
                        if atom_linked.numNitrogenAtoms > 0:
                            if not mol.atoms.elementAt(index).isProtonatedNitrogen:
                                if num_methyl == 3 and mol.atoms.elementAt(index).numHydrogenAtoms == 0:
                                    #  HGAAM0: aliphatic H, NEUTRAL trimethylamine (#)
                                    atom.atomType = "HGAAM0"
                                elif num_methyl == 2 and mol.atoms.elementAt(index).numHydrogenAtoms == 1:
                                    #  HGAAM1: aliphatic H, NEUTRAL dimethylamine (#)
                                    atom.atomType = "HGAAM1"
                                elif num_methyl == 1 and mol.atoms.elementAt(index).numHydrogenAtoms == 2:
                                    #  HGAAM2: aliphatic H, NEUTRAL methylamine (#)
                                    atom.atomType = "HGAAM2"
                        if atom_linked.numHydrogenAtoms == 3 and atom_linked.numNitrogenAtoms == 1:
                            if mol.atoms.elementAt(index).atomType.lower() == "NG3P0".lower():
                                #  HGP5: polar H on quarternary ammonium salt (choline)
                                atom.atomType = "HGP5"
                    if atom_linked.num_linkages == 2 and hasTripleBond(atom_linked):
                        #  HGPAM1: polar H, NEUTRAL dimethylamine (#), terminal alkyne H
                        atom.atomType = "HGPAM1"
                    #  HGR52: Aldehyde H, formamide H (RCOH); nonpolar H, neutral 5-mem planar ring C adjacent to heteroatom or + charge
                    #  CG201, CG202, CG203, CG204
                    #  CG2R51, CG2R52, CG2R53 linked to heteroatom
                    #  CG2D1, CG2DC1 linked to protonated N
                    if atomType != None and (atomType.lower() == "CG2O1".lower() or atomType.lower() == "CG2O2".lower() or atomType.lower() == "CG2O3".lower() or atomType.lower() == "CG2O4".lower()):
                        atom.atomType = "HGR52"
                    elif atomType != None and (atomType.lower() == "CG2R51".lower() or atomType.lower() == "CG2R52".lower() or atomType.lower() == "CG2R53".lower()) and hasHeteroAtom(mol, atom_linked):
                        atom.atomType = "HGR52"
                    elif atomType != None and (atomType.lower() == "CG2D1".lower() or atomType.lower() == "CG2DC1".lower()) and hasProtonatedNitrogen(mol, atom_linked):
                        atom.atomType = "HGR52"
                    if atomType != None and atomType.lower() == "CG2R53".lower() and hasProtonatedNitrogen(mol, atom_linked):
                        #  HGR53: nonpolar H, +ve charge HIS he1(+1)
                        atom.atomType = "HGR53"
                    if atom_linked.isAromatic and atom_linked.numRingAtoms[0] == 5 and atomType.lower() == "CG2R51".lower() and not hasHeteroAtom(mol, atom_linked):
                        #  HGR51: nonpolar H, neutral 5-mem planar ring C, LJ based on benzene
                        atom.atomType = "HGR51"
                    if atom_linked.isRingAtom:
                        if atom_linked.isAromatic:
                            if atom_linked.numRingAtoms[0] == 6:
                                if isTypeHGR62(mol, atom_linked):
                                    #  HGR62: nonpolar H, neutral 6-mem planar ring C adjacent to heteroatom
                                    atom.atomType = "HGR62"
                                if atom_linked.isProtonatedPyridine or atom_linked.isProtonatedPyrimidine:
                                    #  HGR63: nonpolar H, NAD+ nicotineamide all ring CH hydrogens
                                    atom.atomType = "HGR63"
                            elif atom_linked.numRingAtoms[0] == 7:
                                #  HGR71: nonpolar H, neutral 7-mem arom ring, AZUL, azulene, kevo
                                atom.atomType = "HGR71"
                elif element.lower() == "N".lower():
                    if not atom_linked.isAromatic:
                        if not atom_linked.isProtonatedNitrogen:
                            if num_methyl == 2:
                                #  HGAAM1: aliphatic H, NEUTRAL dimethylamine (#)
                                atom.atomType = "HGAAM1"
                            elif num_methyl == 1:
                                #  HGAAM2: aliphatic H, NEUTRAL methylamine (#)
                                atom.atomType = "HGAAM2"
                    if atomType != None and atomType.lower() == "NG311".lower():
                        if atom_linked.numCarbonAtoms == 2:
                            #  HGPAM1: polar H, NEUTRAL dimethylamine (#), terminal alkyne H
                            atom.atomType = "HGPAM1"
                        elif atom_linked.numCarbonAtoms == 1 and (atom_linked.numSulfurAtoms == 1 or atom_linked.numOxygenAtoms == 1):
                            #  HGP1: polar H
                            atom.atomType = "HGP1"
                    elif atomType != None and atomType.lower() == "NG321".lower():
                        if atom_linked.numCarbonAtoms == 1:
                            #  HGPAM2: polar H, NEUTRAL methylamine (#)
                            atom.atomType = "HGPAM2"
                        elif atom_linked.numSulfurAtoms == 1:
                            #  HGP1: polar H
                            atom.atomType = "HGP1"
                    elif atom_linked.num_linkages == 3 and not atom_linked.isProtonatedNitrogen and atom_linked.numHydrogenAtoms == 3:
                        #  HGPAM3: polar H, NEUTRAL ammonia (#)
                        atom.atomType = "HGPAM3"
                    if atomType != None and (atomType.lower() == "NG2P1".lower() or atomType.lower() == "NG3P3".lower() or atomType.lower() == "NG2R52".lower() or atomType.lower() == "NG3P2".lower() or atomType.lower() == "NG3P1".lower()):
                        #  HGP2: polar H, +ve charge
                        #  NG2P1, NG3P3, NG2R52, NG3P2, NG3P1
                        atom.atomType = "HGP2"
                    elif atomType != None and atomType.lower() == "NG2P1".lower():
                        #  HGP3: polar H, thiol
                        atom.atomType = "HGP3"
                elif element.lower() == "S".lower():
                    if atomType != None and atomType.lower() == "SG311".lower():
                        #  HGP3: polar H, thiol
                        atom.atomType = "HGP3"
                if atomType != None and (atomType.lower() == "OG311".lower() or atomType.lower() == "NG2S2".lower() or atomType.lower() == "NG2S1".lower() or atomType.lower() == "NG2R51".lower() or atomType.lower() == "NG2S3".lower() or atomType.lower() == "NG2R61".lower() or atomType.lower() == "NG2D1".lower() or atomType.lower() == "NG3C51".lower() or atomType.lower() == "NG2R53".lower() or atomType.lower() == "NG3N1".lower() or atomType.lower() == "NG2R43".lower()):
                    #  HGP1: polar H
                    #  OG311, NG2S2, NG2S1, NG2R51, NG2S3,
                    #  NG2R61, NG2D1, NG3C51, NG2R53, NG3N1,
                    #  NG2R43
                    atom.atomType = "HGP1"
                    if atomType.lower() == "NG321".lower() and atom_linked.numSulfurAtoms == 1:
                        atom.atomType = "HGP1"
                    if atom_linked.isProtonatedNitrogen and atomType.lower() == "NG2R61".lower():
                        #  HGP2: polar H, +ve charge
                        atom.atomType = "HGP2"
                if atom_linked.element.lower() == "N".lower() and not atom_linked.isProtonatedNitrogen and atom_linked.numHydrogenAtoms == 2 and hasAromaticBond(mol, atom_linked):
                    #  HGP4: polar H, neutral conjugated -NH2 group (NA bases)
                    atom.atomType = "HGP4"
            i += 1

    def hasCarbonDoubleBond(self, mol, atom):
        i = 0
        while i < atom.num_linkages:
            if mol.atoms.elementAt(atom.linkage[i]).element.lower() == "C".lower() and atom.bondOrder[i] == 2:
                return True
            i += 1
        return False

    def hasDoubleBond(self, mol, atm_index):
        atom = mol.atoms.elementAt(atm_index)
        i = 0
        while i < atom.num_linkages:
            if atom.bondOrder[i] == 2:
                return True
            i += 1
        return False

#    def hasDoubleBond(self, atm):
#        i = 0
#        while i < atm.num_linkages:
#            if atm.bondOrder[i] == 2:
#                return True
#            i += 1
#        return False

    def hasTripleBond(self, atm):
        i = 0
        while i < atm.num_linkages:
            if atm.bondOrder[i] == 3:
                return True
            i += 1
        return False

    def countSpecificElement(self, mol, atm, element):
        num_elements = 0
        i = 0
        while i < atm.num_linkages:
            if mol.atoms.elementAt(atm.linkage[i]).element.lower() == element.lower():
                num_elements += 1
            i += 1
        return num_elements

    def getIndex(self, mol, atm, element):
        i = 0
        while i < atm.num_linkages:
            if mol.atoms.elementAt(atm.linkage[i]).element.lower() == element.lower():
                return atm.linkage[i]
            i += 1
        return -1

    def countMethyl(self, mol, atm):
        num_methyl = 0
        i = 0
        while i < atm.num_linkages:
            if mol.atoms.elementAt(atm.linkage[i]).isMethyl:
                num_methyl += 1
            i += 1
        return num_methyl

    def hasAromaticBond(self, mol, atm):
        i = 0
        while i < atm.num_linkages:
            if mol.atoms.elementAt(atm.linkage[i]).isAromatic:
                return True
            i += 1
        return False

    def hasHeteroAtom(self, mol, atm):
        i = 0
        while i < atm.num_linkages:
            if not mol.atoms.elementAt(atm.linkage[i]).element.lower() == "H".lower() and not mol.atoms.elementAt(atm.linkage[i]).element.lower() == "C".lower():
                return True
            i += 1
        return False

    def isTypeHGR62(self, mol, atm):
        atom_linked = Atom()
        atom_linked_linked = Atom()
        hasAdjacentPolarizedGroup = False
        i = 0
        while i < atm.num_linkages:
            atom_linked = mol.atoms.elementAt(atm.linkage[i])
            if atom_linked.element.lower() == "C".lower():
                while j < atom_linked.num_linkages:
                    atom_linked_linked = mol.atoms.elementAt(atom_linked.linkage[j])
                    if not atom_linked_linked.isRingAtom and ((atom_linked_linked.element.lower() == "O".lower() and atom_linked_linked.num_linkages == 1) or isHalogenAtom(atom_linked_linked.element)):
                        hasAdjacentPolarizedGroup = True
                        break
                    j += 1
            i += 1
        if self.hasHeteroAtom(mol, atm) or hasAdjacentPolarizedGroup or atm.ringHasAmide:
            return True
        return False

    def isHalogenAtom(self, element):
        if element.lower() == "CL".lower() or element.lower() == "BR".lower() or element.lower() == "I".lower() or element.lower() == "F".lower():
            return True
        return False

    def hasProtonatedNitrogen(self, mol, atm):
        i = 0
        while i < atm.num_linkages:
            if mol.atoms.elementAt(atm.linkage[i]).isProtonatedNitrogen:
                return True
            i += 1
        return False

