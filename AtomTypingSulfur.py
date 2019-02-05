#!/usr/bin/env python

class AtomTypingSulfur(object):
    def __init__(self):
        """__init__"""

    def setAtomTypeForSulfurs(self, mol):
        num_atoms = len(mol.atoms)
        atom = Atom()
        atom_linked = Atom()
        index = -1
        i = 0
        while i < num_atoms:
            atom = mol.atoms.elementAt(i)
            if atom.element.lower() == "S".lower():
                atom_linked = mol.atoms.elementAt(atom.linkage[0])
                if atom.num_linkages == 2:
                    #  SG311: sulphur, SH, -S-
                    atom.atomType = "SG311"
                if not atom.isRingAtom:
                    if atom.isThiocarbonylS:
                        #  SG2D1: thiocarbonyl S
                        atom.atomType = "SG2D1"
                    elif atom.num_linkages == 2 and atom.numSulfurAtoms == 1 and atom.numCarbonAtoms == 1:
                        #  SG301: sulfur C-S-S-C type
                        atom.atomType = "SG301"
                    elif atom.isDeprotonatedSulfur:
                        #  SG302: thiolate sulfur (-1)
                        atom.atomType = "SG302"
                    elif atom.num_linkages == 3 and atom.numCarbonAtoms == 2 and atom.numOxygenAtoms == 1:
                        #  SG3O3: neutral sulfoxide sulfur
                        atom.atomType = "SG3O3"
                    if atom.numOxygenAtoms >= 2:
                        if charge == -1:
                            #  SG3O1: sulfate -1 sulfur
                            atom.atomType = "SG3O1"
                        elif charge == 0:
                            #  SG3O2: neutral sulfone/sulfonamide sulfur
                            atom.atomType = "SG3O2"
                else:
                    if atom.numOxygenAtoms >= 2 and getCharge(mol, i) == 0:
                        #  SG3O2: neutral sulfone/sulfonamide sulfur
                        atom.atomType = "SG3O2"
                    elif atom.numRingAtoms[0] == 5:
                        #  SG2R50: THIP, thiophene
                        atom.atomType = "SG2R50"
            i += 1

    def getCharge(self, mol, atm_index):
        atom = mol.atoms.elementAt(atm_index)
        atom_linked = Atom()
        num_OG2P1 = 0
        charge = 0
        i = 0
        while i < atom.num_linkages:
            atom_linked = mol.atoms.elementAt(atom.linkage[i])
            if atom_linked.atomType != None and atom_linked.atomType.lower() == "OG2P1".lower():
                num_OG2P1 += 1
            i += 1
        if num_OG2P1 == 3:
            charge = -1
        elif num_OG2P1 == 2:
            charge = 0
        return charge

