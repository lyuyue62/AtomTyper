#!/usr/bin/env python

class AtomTypingHalogen(object):
    def __init__(self):
        """__init__"""
    def setAtomTypeForHalogens(self, mol):
        num_atoms = len(mol.atoms)
        atom = Atom()
        atom_linked = Atom()
        index = -1
        num_element = 0
        i = 0
        while i < num_atoms:
            atom = mol.atoms.elementAt(i)
            if atom.element.lower() == "CL".lower():
                atom_linked = mol.atoms.elementAt(atom.linkage[0])
                if atom_linked.isAromatic and atom_linked.numRingAtoms[0] == 6:
                    #  CLGR1: CHLB, chlorobenzene
                    atom.atomType = "CLGR1"
                else:
                    num_element = countSpecificElement(mol, atom.linkage[0], "CL")
                    if num_element == 1 or num_element == 2:
                        #  CLGA1: CLET, DCLE, chloroethane, 1,1-dichloroethane
                        atom.atomType = "CLGA1"
                    elif num_element == 3:
                        #  CLGA3: TCLE, 1,1,1-trichloroethane
                        atom.atomType = "CLGA3"
            elif atom.element.lower() == "BR".lower():
                atom_linked = mol.atoms.elementAt(atom.linkage[0])
                if atom_linked.isAromatic and atom_linked.numRingAtoms[0] == 6:
                    #  BRGR1: BROB, bromobenzene
                    atom.atomType = "BRGR1"
                else:
                    num_element = countSpecificElement(mol, atom.linkage[0], "BR")
                    if num_element == 1:
                        #  BRGA1: BRET, bromoethane
                        atom.atomType = "BRGA1"
                    elif num_element == 2:
                        #  BRGA2: DBRE, 1,1-dibromoethane
                        atom.atomType = "BRGA2"
                    elif num_element == 3:
                        #  BRGA3: TBRE, 1,1,1-dibromoethane
                        atom.atomType = "BRGA3"
            elif atom.element.lower() == "I".lower():
                atom_linked = mol.atoms.elementAt(atom.linkage[0])
                if atom_linked.isAromatic and atom_linked.numRingAtoms[0] == 6:
                    #  IGR1: IODB, iodobenzene
                    atom.atomType = "IGR1"
            elif atom.element.lower() == "F".lower():
                atom_linked = mol.atoms.elementAt(atom.linkage[0])
                if atom_linked.isAromatic:
                    #  FGR1: aromatic flourine
                    atom.atomType = "FGR1"
                else:
                    num_element = countSpecificElement(mol, atom.linkage[0], "F")
                    if num_element == 1:
                        #  FGA1: aliphatic fluorine, monofluoro
                        atom.atomType = "FGA1"
                    elif num_element == 2:
                        #  FGA2: aliphatic fluorine, difluoro
                        atom.atomType = "FGA2"
                    elif num_element == 3:
                        #  FGA3: aliphatic fluorine, trifluoro
                        atom.atomType = "FGA3"
                    elif num_element == 4:
                        #  FGP1: anionic F, for ALF4 AlF4-
                        atom.atomType = "FGP1"
            i += 1

    def countSpecificElement(self, mol, atm_index, element):
        atom = mol.atoms.elementAt(atm_index)
        num_elements = 0
        i = 0
        while i < atom.num_linkages:
            if mol.atoms.elementAt(atom.linkage[i]).element.lower() == element.lower():
                num_elements += 1
            i += 1
        return num_elements

