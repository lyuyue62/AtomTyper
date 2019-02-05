#!/usr/bin/env python

class AtomTypingCarbon(object):
    def __init__(self):
        """__init__"""

    def setAtomTypeForCarbons(self, mol, aRingAromacity, aRingHasCabonyl):
        num_atoms = len(mol.atoms)
        atom = Atom()
        i = 0
        while i < num_atoms:
            atom = mol.atoms.elementAt(i)
            if atom.element.lower() == "C".lower():
                if atom.num_linkages == 2:
                    if isInternalAlkyne(mol, i):
                        #  CG1T1: internal alkyne R-C#C
                        atom.atomType = "CG1T1"
                    elif isTerminalAlkyne(mol, i):
                        #  CG1T2: terminal alkyne H-C#C
                        atom.atomType = "CG1T2"
                    elif isForCyanoGroup(mol, i):
                        #  CG1N1: C for cyano group
                        atom.atomType = "CG1N1"
                    elif isCO2Carbon(mol, i):
                        #  CG2O7: CO2 carbon
                        atom.atomType = "CG2O7"
                elif atom.num_linkages == 3:
                    if not atom.isConjugated:
                        if isInternalAlkene(mol, i) or isImine(mol, i):
                            #  CG2D1: internal alkene; RHC= ; imine C (-N=C-)
                            atom.atomType = "CG2D1"
                        elif isTerminalAlkene(mol, i):
                            #  CG2D2: terminal alkene; H2C=
                            atom.atomType = "CG2D2"
                        #  CG2D1O: double bond C adjacent to heteroatom. In conjugated systems, the atom to which it is double bonded must be CG2DC1.
                        if isDoubleBondCAdjacentToHeteroatom(mol, i):
                            atom.atomType = "CG2D1O"
                        # *CG2D2O: double bond C adjacent to heteroatom. In conjugated systems, the atom to which it is double bonded must be CG2DC2.
                        if atom.numNitrogenAtoms == 3:
                            #  CG2N1 : conjugated C in guanidine/guanidinium (carbon linked to three N)
                            atom.atomType = "CG2N1"
                        elif atom.numNitrogenAtoms == 2:
                            #  CG2N2 : conjugated C in amidinium cation (carbon linked to two N)
                            atom.atomType = "CG2N2"
                    else:
                        if isInternalAlkene(mol, i):
                            #  CG2DC1: conjugated alkenes, R2C=CR2
                            atom.atomType = "CG2DC1"
                            # *CG2DC2: conjugated alkenes, R2C=CR2. What is the different between CG2DC1 and CG2DC2?
                        elif atom.numHydrogenAtoms == 2:
                            #  CG2DC3: conjugated alkenes, H2C=
                            atom.atomType = "CG2DC3"
                        #  CG2D1O: double bond C adjacent to heteroatom. In conjugated systems, the atom to which it is double bonded must be CG2DC1.
                        if isDoubleBondCAdjacentToHeteroatom(mol, i):
                            atom.atomType = "CG2D1O"
                        # *CG2D2O: double bond C adjacent to heteroatom. In conjugated systems, the atom to which it is double bonded must be CG2DC2.
                        if atom.numNitrogenAtoms == 3:
                            #  CG2N1 : conjugated C in guanidine/guanidinium (carbon linked to three N)
                            atom.atomType = "CG2N1"
                        elif atom.numNitrogenAtoms == 2:
                            #  CG2N2 : conjugated C in amidinium cation (carbon linked to two N)
                            atom.atomType = "CG2N2"
                    if atom.isCarbonyl:
                        if atom.isAmide:
                            #  CG2O1: carbonyl C: amides (O=C-N)
                            atom.atomType = "CG2O1"
                        elif atom.numOxygenAtoms == 2 and atom.numNitrogenAtoms != 1:
                            if isNeutralCarboxylicAcid(mol, i):
                                #  CG2O2: carbonyl C: esters, [neutral] carboxylic acids
                                atom.atomType = "CG2O2"
                            else:
                                #  CG2O3: carbonyl C: [negative] carboxylates
                                atom.atomType = "CG2O3"
                        elif atom.numOxygenAtoms == 1 and atom.numHydrogenAtoms == 1:
                            #  CG2O4: carbonyl C: aldehydes
                            atom.atomType = "CG2O4"
                        elif atom.isKetone:
                            #  CG2O5: carbonyl C: ketones
                            atom.atomType = "CG2O5"
                        elif isUreaOrCarbonate(mol, i):
                            #  CG2O6: carbonyl C: urea, carbonate (CS3, CNO2, CO3)
                            atom.atomType = "CG2O6"
                    if atom.isRingAtom:
                        if not is_ring_brdiging_atom:
                            if atom.numRingAtoms[0] == 5:
                                if atom.isAromatic:
                                    if atom.atomType != None and atom.atomType.lower() == "CG2DC1".lower() and hasExocyclicDoubleBond(mol, i):
                                        #  CG25C1: same as CG2DC1 but in 5-membered ring with exocyclic double bond
                                        atom.atomType = "CG25C1"
                                        # *CG25C2: same as CG2DC2 but in 5-membered ring with exocyclic double bond
                                    elif atom.atomType != None and atom.atomType.lower() == "CG2D1O".lower() and hasExocyclicDoubleBond(mol, i):
                                        #  CG251O: same as CG2D1O but in 5-membered ring with exocyclic double bond
                                        atom.atomType = "CG251O"
                                        # *CG252O: same as CG2D2O but in 5-membered ring with exocyclic double bond
                                    else:
                                        #  CG2R51: 5-mem ring, his CG, CD2(0), trp
                                        atom.atomType = "CG2R51"
                                        if hasDoubleBondedNitrogen(mol, i):
                                            #  CG2R52: 5-mem ring, double bonded to N, PYRZ, pyrazole
                                            atom.atomType = "CG2R52"
                                        if isTypeCG2R53(mol, i):
                                            #  CG2R53: 5-mem ring, double bonded to heteroatom and adjacent to another heteroatom, purine C8, his CE1 (0,+1), 2PDO, kevo
                                            atom.atomType = "CG2R53"
                                        if isBipyrroleCarbon(mol, i):
                                            #  CG2R57: 5-mem ring, bipyrroles
                                            atom.atomType = "CG2R57"
                                else:
                                    #  CG2R51: 5-mem ring, his CG, CD2(0), trp
                                    atom.atomType = "CG2R51"
                                    if isTypeCG2R53(mol, i):
                                        #  CG2R53: 5-mem ring, double bonded to heteroatom and adjacent to another heteroatom, purine C8, his CE1 (0,+1), 2PDO, kevo
                                        atom.atomType = "CG2R53"
                            elif atom.numRingAtoms[0] == 6:
                                if atom.isAromatic:
                                    #  CG2R61: 6-mem aromatic C
                                    atom.atomType = "CG2R61"
                                    if atom.isProtonatedPyridine or atom.isProtonatedPyrimidine or (atom.ringHasCabonyl and not atom.isCarbonyl):
                                        #  CG2R62: 6-mem aromatic C for protonated pyridine (NIC) and rings containing carbonyls (see CG2R63) (NA)
                                        atom.atomType = "CG2R62"
                                    if atom.isCarbonyl:
                                        #  CG2R63: 6-mem aromatic amide carbon (NA) (and other 6-mem aromatic carbonyls?)
                                        atom.atomType = "CG2R63"
                                    if not atom.isCarbonyl and isTypeCG2R64(mol, i):
                                        #  CG2R64: 6-mem aromatic amidine and guanidine carbon (between 2 or 3 Ns and double-bound to one of them), NA, PYRM
                                        atom.atomType = "CG2R64"
                                    if isTypeCG2R66(mol, i):
                                        #  CG2R66: 6-mem aromatic carbon bound to F
                                        atom.atomType = "CG2R66"
                                    if isTypeCG2R67(mol, i, is_ring_brdiging_atom, aRingAromacity):
                                        #  CG2R67: 6-mem aromatic carbon of biphenyl
                                        atom.atomType = "CG2R67"
                            elif atom.numRingAtoms[0] == 7:
                                if atom.isAromatic:
                                    #  CG2R71: 7-mem ring arom C, AZUL, azulene, kevo
                                    atom.atomType = "CG2R71"
                        else:
                            if ring_index_5 != -1:
                                #  CG2R51: 5-mem ring, his CG, CD2(0), trp
                                atom.atomType = "CG2R51"
                            if ring_index_5 != -1 and isTypeCG2R53(mol, i):
                                #  CG2R53: 5-mem ring, double bonded to heteroatom and adjacent to another heteroatom, purine C8, his CE1 (0,+1), 2PDO, kevo
                                atom.atomType = "CG2R53"
                            if ring_index_6 != -1 and bridgedAtomBelongsToAromaticRing(mol, i, 6, aRingAromacity):
                                #  CG2R61: 6-mem aromatic C
                                atom.atomType = "CG2R61"
                                if atom.isProtonatedPyridine or atom.isProtonatedPyrimidine:
                                    #  CG2R62: 6-mem aromatic C for protonated pyridine (NIC) and rings containing carbonyls (see CG2R63) (NA)
                                    atom.atomType = "CG2R62"
                                elif atom.ringHasCabonyl and not atom.isCarbonyl:
                                    while j < 3:
                                        if atom.numRingAtoms[j] == 6 and aRingAromacity[atom.ringIndex[j]] and aRingHasCabonyl[atom.ringIndex[j]]:
                                            #  CG2R62: 6-mem aromatic C for protonated pyridine (NIC) and rings containing carbonyls (see CG2R63) (NA)
                                            atom.atomType = "CG2R62"
                                        j += 1
                                if not atom.isCarbonyl and isTypeCG2R64(mol, i):
                                    #  CG2R64: 6-mem aromatic amidine and guanidine carbon (between 2 or 3 Ns and double-bound to one of them), NA, PYRM
                                    atom.atomType = "CG2R64"
                                if isTypeCG2R67(mol, i, is_ring_brdiging_atom, aRingAromacity):
                                    #  CG2R67: 6-mem aromatic carbon of biphenyl
                                    atom.atomType = "CG2R67"
                            if isTypeCG2RC0(mol, i, aRingAromacity):
                                if atom.atomType == None or not atom.atomType.lower() == "CG2R67".lower():
                                    #  CG2RC0: 6/5-mem ring bridging C, guanine C4,C5, trp
                                    atom.atomType = "CG2RC0"
                            if ring_index_7 != -1 and bridgedAtomBelongsToAromaticRing(mol, i, 7, aRingAromacity):
                                #  CG2RC7: sp2 ring connection with single bond(!), AZUL, azulene, kevo
                                atom.atomType = "CG2RC7"
                elif atom.num_linkages == 4:
                    if atom.numHydrogenAtoms == 0:
                        #  CG301: aliphatic C, no hydrogens, neopentane
                        atom.atomType = "CG301"
                    elif atom.numHydrogenAtoms == 1:
                        #  CG311: aliphatic C with 1 H, CH
                        atom.atomType = "CG311"
                        if hasAdjecentProtonatedNitrogen(mol, i):
                            #  CG314: aliphatic C with 1 H, adjacent to positive N (PROT NTER) (+)
                            atom.atomType = "CG314"
                    elif atom.numHydrogenAtoms == 2:
                        #  CG321: aliphatic C for CH2
                        atom.atomType = "CG321"
                        if hasAdjecentProtonatedNitrogen(mol, i):
                            #  CG324: aliphatic C for CH2, adjacent to positive N (piperidine) (+)
                            atom.atomType = "CG324"
                    elif atom.numHydrogenAtoms == 3:
                        #  CG331: aliphatic C for methyl group (-CH3)
                        atom.atomType = "CG331"
                        if index != -1:
                            if mol.atoms.elementAt(index).isProtonatedNitrogen:
                                #  CG334: aliphatic C for methyl group (-CH3), adjacent to positive N (PROT NTER) (+)
                                atom.atomType = "CG334"
                            else:
                                if not hasDoubleBond(mol, index):
                                    if countMethyl(mol, index) == 3:
                                        #  CG3AM0: aliphatic C for CH3, NEUTRAL trimethylamine methyl carbon (#)
                                        atom.atomType = "CG3AM0"
                                    elif countMethyl(mol, index) == 2 and mol.atoms.elementAt(index).numHydrogenAtoms == 1:
                                        #  CG3AM1: aliphatic C for CH3, NEUTRAL dimethylamine methyl carbon (#)
                                        atom.atomType = "CG3AM1"
                                    elif countMethyl(mol, index) == 1 and mol.atoms.elementAt(index).numHydrogenAtoms == 2:
                                        #  CG3AM2: aliphatic C for CH3, NEUTRAL methylamine methyl carbon (#)
                                        atom.atomType = "CG3AM2"
                    if num_F == 3:
                        #  CG302: aliphatic C, no hydrogens, trifluoromethyl
                        atom.atomType = "CG302"
                    elif num_F == 2:
                        #  CG312: aliphatic C with 1 H, difluoromethyl
                        atom.atomType = "CG312"
                    elif num_F == 1:
                        #  CG322: aliphatic C for CH2, monofluoromethyl
                        atom.atomType = "CG322"
                    if atom.numSulfurAtoms == 1:
                        #  CG323: aliphatic C for CH2, thiolate carbon
                        if index != -1 and mol.atoms.elementAt(index).isDeprotonatedSulfur:
                            atom.atomType = "CG323"
                    if atom.isRingAtom:
                        if not is_ring_brdiging_atom:
                            if atom.numRingAtoms[0] == 3:
                                #  CG3C31: cyclopropyl carbon
                                atom.atomType = "CG3C31"
                            elif atom.numRingAtoms[0] == 4:
                                #  CG3C41: cyclobutyl carbon
                                atom.atomType = "CG3C41"
                            elif atom.numRingAtoms[0] == 5:
                                if atom.numHydrogenAtoms == 0:
                                    #  CG3C50: 5-mem ring aliphatic quaternary C (cholesterol, bile acids)
                                    atom.atomType = "CG3C50"
                                elif atom.numHydrogenAtoms == 1:
                                    #  CG3C51: 5-mem ring aliphatic CH  (proline CA, furanoses)
                                    atom.atomType = "CG3C51"
                                    #  CG3C53: 5-mem ring aliphatic CH  adjacent to positive N (proline.H+ CA) (+)
                                    if hasAdjecentProtonatedNitrogen(mol, i):
                                        atom.atomType = "CG3C53"
                                elif atom.numHydrogenAtoms == 2:
                                    #  CG3C52: 5-mem ring aliphatic CH2 (proline CB/CG/CD, THF, deoxyribose)
                                    atom.atomType = "CG3C52"
                                    #  CG3C54: 5-mem ring aliphatic CH2 adjacent to positive N (proline.H+ CD) (+)
                                    if hasAdjecentProtonatedNitrogen(mol, i):
                                        atom.atomType = "CG3C54"
                        else:
                            if getAtomsOfSmallestRing(mol, i) <= 5:
                                #  CG3RC1: bridgehead in bicyclic systems containing at least one 5-membered or smaller ring
                                atom.atomType = "CG3RC1"
            i += 1

    def getIndex(self, mol, atm_index, element):
        atom = mol.atoms.elementAt(atm_index)
        i = 0
        while i < atom.num_linkages:
            if mol.atoms.elementAt(atom.linkage[i]).element.lower() == element.lower():
                return atom.linkage[i]
            i += 1
        return -1

    def countSpecificElement(self, mol, atm_index, element):
        atom = mol.atoms.elementAt(atm_index)
        num_elements = 0
        i = 0
        while i < atom.num_linkages:
            if mol.atoms.elementAt(atom.linkage[i]).element.lower() == element.lower():
                num_elements += 1
            i += 1
        return num_elements

    def hasDoubleBond(self, mol, atm_index):
        atom = mol.atoms.elementAt(atm_index)
        i = 0
        while i < atom.num_linkages:
            if atom.bondOrder[i] == 2:
                return True
            i += 1
        return False

    def countMethyl(self, mol, atm_index):
        atom = mol.atoms.elementAt(atm_index)
        num_methyls = 0
        i = 0
        while i < atom.num_linkages:
            if mol.atoms.elementAt(atom.linkage[i]).isMethyl:
                num_methyls += 1
            i += 1
        return num_methyls

    def getAtomsOfSmallestRing(self, mol, atm_index):
        atom = mol.atoms.elementAt(atm_index)
        min_num_ring_atoms = 100
        i = 0
        while i < 3:
            if atom.ringIndex[i] != -1 and atom.numRingAtoms[i] < min_num_ring_atoms:
                min_num_ring_atoms = atom.numRingAtoms[i]
            i += 1
        return min_num_ring_atoms

    def isAtomInFiveMemberedRing(self, mol, atm_index):
        atom = mol.atoms.elementAt(atm_index)
        i = 0
        while i < 3:
            if atom.numRingAtoms[i] == 5:
                return atom.ringIndex[i]
            i += 1
        return -1

    def isAtomInSixMemberedRing(self, mol, atm_index):
        atom = mol.atoms.elementAt(atm_index)
        i = 0
        while i < 3:
            if atom.numRingAtoms[i] == 6:
                return atom.ringIndex[i]
            i += 1
        return -1

    def isAtomInSevenMemberedRing(self, mol, atm_index):
        atom = mol.atoms.elementAt(atm_index)
        i = 0
        while i < 3:
            if atom.numRingAtoms[i] == 7:
                return atom.ringIndex[i]
            i += 1
        return -1

    def hasDifferentRingIndex(self, mol, atom_i, atom_j):
        index_i = int()
        index_j = int()
        i = 0
        while i < 3:
            index_i = mol.atoms.elementAt(atom_i).ringIndex[i]
            while j < 3:
                index_j = mol.atoms.elementAt(atom_j).ringIndex[j]
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
            index_i = mol.atoms.elementAt(atom_i).ringIndex[i]
            while j < 3:
                index_j = mol.atoms.elementAt(atom_j).ringIndex[j]
                if index_i != -1 and index_j != -1 and index_i == index_j:
                    return True
                j += 1
            i += 1
        return False

    def hasExocyclicDoubleBond(self, mol, atm_index):
        atom = mol.atoms.elementAt(atm_index)
        atom_linked = Atom()
        if atom.isRingAtom:
            while i < atom.num_linkages:
                atom_linked = mol.atoms.elementAt(atom.linkage[i])
                if not self.hasSameRingIndex(mol, atm_index, atom.linkage[i]) and atom.bondOrder[i] == 2:
                    return True
                i += 1
        return False

    def bridgedAtomBelongsToAromaticRing(self, mol, atm_index, num_ring_atoms, aRingAromacity):
        atom = mol.atoms.elementAt(atm_index)
        i = 0
        while i < 3:
            if atom.ringIndex[i] != -1 and atom.numRingAtoms[i] == num_ring_atoms and aRingAromacity[atom.ringIndex[i]]:
                return True
            i += 1
        return False

    def allBridgedRingsAreAromatic(self, mol, atm_index, aRingAromacity):
        atom = mol.atoms.elementAt(atm_index)
        i = 0
        while i < 3:
            if atom.ringIndex[i] != -1 and not aRingAromacity[atom.ringIndex[i]]:
                return False
            i += 1
        return True

    def isInternalAlkyne(self, mol, atm_index):
        atom = mol.atoms.elementAt(atm_index)
        i = 0
        while i < atom.num_linkages:
            if atom.bondOrder[i] == 3 and mol.atoms.elementAt(atom.linkage[i]).element.lower() == "C".lower():
                if atom.numHydrogenAtoms == 0:
                    return True
            i += 1
        return False

    def isTerminalAlkyne(self, mol, atm_index):
        atom = mol.atoms.elementAt(atm_index)
        i = 0
        while i < atom.num_linkages:
            if atom.bondOrder[i] == 3 and mol.atoms.elementAt(atom.linkage[i]).element.lower() == "C".lower():
                if atom.numHydrogenAtoms == 1:
                    return True
            i += 1
        return False

    def isForCyanoGroup(self, mol, atm_index):
        atom = mol.atoms.elementAt(atm_index)
        i = 0
        while i < atom.num_linkages:
            if atom.bondOrder[i] == 3 and mol.atoms.elementAt(atom.linkage[i]).element.lower() == "N".lower():
                return True
            i += 1
        return False

    def isInternalAlkene(self, mol, atm_index):
        atom = mol.atoms.elementAt(atm_index)
        i = 0
        while i < atom.num_linkages:
            if atom.bondOrder[i] == 2 and (mol.atoms.elementAt(atom.linkage[i]).element.lower() == "C".lower() or mol.atoms.elementAt(atom.linkage[i]).element.lower() == "N".lower()):
                if atom.numHydrogenAtoms != 2:
                    return True
            i += 1
        return False

    def isDoubleBondCAdjacentToHeteroatom(self, mol, atm_index):
        atom = mol.atoms.elementAt(atm_index)
        i = 0
        while i < atom.num_linkages:
            if atom.bondOrder[i] == 2 and mol.atoms.elementAt(atom.linkage[i]).element.lower() == "C".lower():
                while j < atom.num_linkages:
                    if j != i:
                        if not mol.atoms.elementAt(atom.linkage[j]).element.lower() == "H".lower() and not mol.atoms.elementAt(atom.linkage[j]).element.lower() == "C".lower():
                            return True
                    j += 1
            i += 1
        return False

    def isImine(self, mol, atm_index):
        atom = mol.atoms.elementAt(atm_index)
        i = 0
        while i < atom.num_linkages:
            if atom.bondOrder[i] == 2 and mol.atoms.elementAt(atom.linkage[i]).element.lower() == "N".lower():
                return True
            i += 1
        return False

    def isTerminalAlkene(self, mol, atm_index):
        atom = mol.atoms.elementAt(atm_index)
        i = 0
        while i < atom.num_linkages:
            if atom.bondOrder[i] == 2 and mol.atoms.elementAt(atom.linkage[i]).element.lower() == "C".lower():
                if atom.numHydrogenAtoms == 2:
                    return True
            i += 1
        return False

    def isNeutralCarboxylicAcid(self, mol, atm_index):
        atom = mol.atoms.elementAt(atm_index)
        atom_linked = Atom()
        if atom.numOxygenAtoms == 2:
            while i < atom.num_linkages:
                atom_linked = mol.atoms.elementAt(atom.linkage[i])
                if atom.bondOrder[i] == 1 and atom_linked.element.lower() == "O".lower() and atom_linked.num_linkages == 2:
                    return True
                i += 1
        return False

    def isUreaOrCarbonate(self, mol, atm_index):
        atom = mol.atoms.elementAt(atm_index)
        if atom.numSulfurAtoms == 3:
            return True
        if atom.numOxygenAtoms == 3:
            return True
        if atom.numNitrogenAtoms == 1 and atom.numOxygenAtoms == 2:
            return True
        return False

    def isCO2Carbon(self, mol, atm_index):
        atom = mol.atoms.elementAt(atm_index)
        if atom.num_linkages == 2 and (atom.numOxygenAtoms == 2 or (atom.numOxygenAtoms == 1 and atom.numNitrogenAtoms == 1)):
            return True
        return False

    def hasDoubleBondedNitrogen(self, mol, atm_index):
        atom = mol.atoms.elementAt(atm_index)
        i = 0
        while i < 3:
            if atom.linkage[i] != -1 and mol.atoms.elementAt(atom.linkage[i]).element.lower() == "N".lower() and atom.bondOrder[i] == 2:
                return True
            i += 1
        return False

    def isTypeCG2R53(self, mol, atm_index):
        atom = mol.atoms.elementAt(atm_index)
        atom_linked = Atom()
        hasHeteroAtomSingleBond = False
        hasHeteroAtomDoubleBond = False
        i = 0
        while i < atom.num_linkages:
            atom_linked = mol.atoms.elementAt(atom.linkage[i])
            if atom.bondOrder[i] == 1 and not atom_linked.element.lower() == "H".lower() and not atom_linked.element.lower() == "C".lower():
                hasHeteroAtomSingleBond = True
            if atom.bondOrder[i] == 2 and not atom_linked.element.lower() == "H".lower() and not atom_linked.element.lower() == "C".lower():
                hasHeteroAtomDoubleBond = True
            i += 1
        if hasHeteroAtomSingleBond and hasHeteroAtomDoubleBond:
            return True
        return False

    def isBipyrroleCarbon(self, mol, atm_index):
        atom = mol.atoms.elementAt(atm_index)
        atom_linked = Atom()
        if atom.isPyrrole:
            while i < atom.num_linkages:
                atom_linked = mol.atoms.elementAt(atom.linkage[i])
                if atom_linked.isPyrrole and not mol.atoms.elementAt(atom.linkage[i]).isBridgingAtom and (atom.ringIndex[0] != atom_linked.ringIndex[0]):
                    return True
                i += 1
        return False

    def isTypeCG25C1(self, mol, atm_index):
        atom = mol.atoms.elementAt(atm_index)
        i = 0
        while i < atom.num_linkages:
            if atom.bondOrder[i] == 2 and mol.atoms.elementAt(atom.linkage[i]).element.lower() == "C".lower():
                if atom.numHydrogenAtoms != 2:
                    return True
            i += 1
        return False

    def isTypeCG2R64(self, mol, atm_index):
        # between 2 or 3 Ns and double-bound to one of them
        atom = mol.atoms.elementAt(atm_index)
        if atom.numNitrogenAtoms == 2 or atom.numNitrogenAtoms == 3:
            while i < atom.num_linkages:
                if mol.atoms.elementAt(atom.linkage[i]).element.lower() == "N".lower() and mol.atoms.elementAt(atom.linkage[i]).isRingAtom and (atom.bondOrder[i] == 2 or mol.atoms.elementAt(atom.linkage[i]).isAromatic):
                    return True
                i += 1
        return False

    def isTypeCG2R66(self, mol, atm_index):
        atom = mol.atoms.elementAt(atm_index)
        i = 0
        while i < 4:
            if atom.linkage[i] != -1 and mol.atoms.elementAt(atom.linkage[i]).element.lower() == "F".lower():
                return True
            i += 1
        return False

    def isTypeCG2R67(self, mol, atm_index, is_ring_bridging_atom, aRingAromacity):
        atom = mol.atoms.elementAt(atm_index)
        atom_linked = Atom()
        index_i = -1
        index_j = -1
        if not is_ring_bridging_atom:
            if atom.numRingAtoms[0] == 6 and atom.isAromatic:
                while i < atom.num_linkages:
                    atom_linked = mol.atoms.elementAt(atom.linkage[i])
                    if not mol.atoms.elementAt(atom.linkage[i]).isBridgingAtom and atom_linked.numRingAtoms[0] == 6 and atom_linked.isAromatic:
                        if numOverlappedRings(mol, atom, atom_linked) == 0:
                            return True
                    i += 1
        else:
            if ringIndex != -1 and aRingAromacity[ringIndex]:
                while i < atom.num_linkages:
                    atom_linked = mol.atoms.elementAt(atom.linkage[i])
                    ringIndex = self.isAtomInSixMemberedRing(mol, atom.linkage[i])
                    if mol.atoms.elementAt(atom.linkage[i]).isBridgingAtom and ringIndex != -1 and aRingAromacity[ringIndex]:
                        if numOverlappedRings(mol, atom, atom_linked) == 1:
                            return True
                    i += 1
        return False

    def numOverlappedRings(self, mol, atm1, atm2):
        num_overlapped_rings = 0
        i = 0
        while i < 3:
            while j < 3:
                if atm1.ringIndex[i] != -1 and atm2.ringIndex[j] != -1 and atm1.ringIndex[i] == atm2.ringIndex[j]:
                    num_overlapped_rings += 1
                j += 1
            i += 1
        return num_overlapped_rings

    def isTypeCG2RC0(self, mol, atm_index, aRingAromacity):
        atom = mol.atoms.elementAt(atm_index)
        fiveMemberedRingIndex = self.isAtomInFiveMemberedRing(mol, atm_index)
        sixMemberedRingIndex = self.isAtomInSixMemberedRing(mol, atm_index)
        if fiveMemberedRingIndex != -1 and sixMemberedRingIndex != -1 and aRingAromacity[sixMemberedRingIndex]:
            return True
        return False

    def hasAdjecentProtonatedNitrogen(self, mol, atm_index):
        atom = mol.atoms.elementAt(atm_index)
        atom_linked = Atom()
        atom_linked_linked = Atom()
        i = 0
        while i < atom.num_linkages:
            atom_linked = mol.atoms.elementAt(atom.linkage[i])
            if atom_linked.element.lower() == "N".lower():
                if atom_linked.isProtonatedNitrogen:
                    return True
                else:
                    while j < atom_linked.num_linkages:
                        if atom_linked.linkage[j] != atm_index:
                            atom_linked_linked = mol.atoms.elementAt(atom_linked.linkage[j])
                            if atom_linked_linked.isProtonatedImineGroup:
                                return True
                        j += 1
            i += 1
        return False

