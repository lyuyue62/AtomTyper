#!/usr/bin/env python

import math
import traceback
import logging

from Edge import Edge
from Atom import Atom
from Bond import Bond
from Fragment import Fragment
from SmallMolecule import SmallMolecule
from LinkedAtomGroup import LinkedAtomGroup


class ChemicalToolKits(object):
    def __init__(self):
        """ __init__ """
    def getDistance(self, r1, r2):
        return math.sqrt((r1[0] - r2[0]) * (r1[0] - r2[0]) + (r1[1] - r2[1]) * (r1[1] - r2[1]) + (r1[2] - r2[2]) * (r1[2] - r2[2]))

    def readMol2File(self, file_name):
        mol = SmallMolecule()
        atom = Atom()
        bond = Bond()
        line = str()
        element = str()
        index = int()
        i = int()
        j = int()

        with open (file_name) as f:
            for line in f:
                #  *** read atom information ***
                if line.startswith('@<TRIPOS>ATOM'):
                    for line in f:
                        if line.startswith('@<TRIPOS>BOND'):
                            break
                        #  atom_id / atom_name / x / y / z / atom_type / subst_id / subst_name / charge / status_bit
                        atom = Atom()
                        atom.line = line
                        atom.atom_name = line.strip().split()[1]
                        atom.R[0] = float(line.strip().split()[2])
                        atom.R[1] = float(line.strip().split()[3])
                        atom.R[2] = float(line.strip().split()[4])
                        element = line.strip().split()[5]
                        index = element.find('.')
                        if index != -1:
                            atom.element = element[0: index]
                        else:
                            atom.element = element
                        atom.element = atom.element.upper()
                        mol.atoms.append(atom)

                #  *** read bond information ***
                if line.startswith('@<TRIPOS>BOND'):
                    for line in f:
                        if line.startswith('@<TRIPOS>SUBSTRUCTURE') or (0 == len(line)):
                            break
                        #  bond_id / origin_atom_id / target_atom_id / bond_type / status_bits
                        #  Bond Types:
                        #  1 = single
                        #  2 = double
                        #  3 = triple
                        #  am = amide
                        #  ar = aromatic
                        #  du = dummy
                        #  un = unknown (cannot be determined from the parameter tables)
                        #  nc = not connected
                        bond = Bond()
                        i = int(line.strip().split()[1]) - 1
                        j = int(line.strip().split()[2]) - 1
                        if i < j:
                            bond.i = i
                            bond.j = j
                        else:
                            bond.i = j
                            bond.j = i
                        bond.bondType = line.strip().split()[3]
                        mol.bonds.append(bond)
                    break

        f.close()
        return mol

    def readCGenFFBondParameters(self, file_name):
        hBondParameter = {}
        try:
            with open (file_name) as f:
                for line in f:
                    if line.startswith('BONDS'):
                        for line in f:
                            if len(line) == 0 or line.startswith('ANGLES'):
                                break
                            if line[0] != '!':
                                st = line.split()
                                hBondParameter[st[0] + " " + st[1]] = st[2] + " " + st[3]
                    break
            f.close()

        except Exception as e:
            logging.error(traceback.format_exc())
        return hBondParameter

    def readCGenFFAngleParameters(self, file_name):
        hAngleParameter = {}
        try:
            with open (file_name) as f:
                for line in f:
                    if line.startswith('ANGLES'):
                        for line in f:
                            if len(line) == 0 or line.startswith('DIHEDRALS'):
                                break
                            if line[0] != '!':
                                st = line.split()
                                hBondParameter[st[0] + " " + st[1] + " " + st[2]] = st[3] + " " + st[4]
                    break
            f.close()

        except Exception as e:
            logging.error(traceback.format_exc())
        return hAngleParameter

    def readCGenFFDihedralParameters(self, file_name):
        hDihedralParameter = {}
        try:
            with open (file_name) as f:
                for line in f:
                    if line.startswith('DIHEDRALS'):
                        for line in f:
                            if len(line) == 0 or line.startswith('IMPROPERS'):
                                break
                            if line[0] != '!':
                                st = line.split()
                                hBondParameter[st[0] + " " + st[1] + " " + st[2] + " " + st[3]] = st[4] + " " + st[5] + " " + st[6]
                    break
            f.close()

        except Exception as e:
            logging.error(traceback.format_exc())
        return hDihedralParameter

    def readCGenFFImproperParameters(self, file_name):
        hImproperParameter = {}
        try:
            with open (file_name) as f:
                for line in f:
                    if line.startswith('IMPROPERS'):
                        for line in f:
                            if len(line) == 0 or line.startswith('NONBONDED'):
                                break
                            if line[0] != '!':
                                st = line.split()
                                hBondParameter[st[0] + " " + st[1] + " " + st[2] + " " + st[3]] = st[4] + " " + st[5] + " " + st[6]
                    break
            f.close()

        except Exception as e:
            logging.error(traceback.format_exc())
        return hImproperParameter

    def readAtomPriority(self, file_name):
        hAtomPriority = {}
        try:
            with open (file_name) as f:
                for line in f:
                    array = line.strip().split("\t")
                    hAtomPriority[array[0]] = array[1]
            f.close()

        except Exception as e:
            logging.error(traceback.format_exc())
        return hAtomPriority

    def readAtomTypes(self, file_name):
        hAtomType = {}
        try:
            with open (file_name) as f:
                for line in f:
                    hAtomType[line[10:16].strip()] = line[29, line.find('!')].strip()
            f.close()

        except Exception as e:
            logging.error(traceback.format_exc())
        return hAtomType

    def setElement(self, mol, hAtomType):
        num_atoms = len(mol.atoms)
        for i in range(num_atoms):
            mol.atoms[i].element = hAtomType[mol.atoms[i].atomType]

    def readChargeIncrementBond(self, file_name):
        hChargeIncrementBond = {}
        try:
            with open (file_name) as f:
                for line in f:
                    array = line.split("\t")
                    hChargeIncrementBond[array[0]] = array[1]
            f.close()

        except Exception as e:
            logging.error(traceback.format_exc())
        return hChargeIncrementBond

    def orderAtoms(self, hAtomPriority, atoms_combi, kind):
        atom_priority_i = 0
        atom_priority_j = 0
        atom_priority_k = 0
        atom_priority_l = 0

        if kind == 1:
        # BONDS: THE LOWEST PRIORITY ATOM ALWAYS COMES FIRST
            array = atoms_combi.split(" ")
            atom_priority_i = int(hAtomPriority[array[0]])
            atom_priority_j = int(hAtomPriority[array[1]])
            if atom_priority_i <= atom_priority_j:
                return atoms_combi
            else:
                return array[1] + " " + array[0]

        elif kind == 2:
        # ANGLES: IF COLUMN 3 HAS A LOWER PRIORITY THAN COLUMN 1, COLUMNS 1 & 3 ARE SWAPPED
            array = atoms_combi.split(" ")
            atom_priority_i = int(hAtomPriority[array[0]])
            atom_priority_k = int(hAtomPriority[array[2]])
            if atom_priority_i <= atom_priority_k:
                return atoms_combi
            else:
                return array[2] + " " + array[1] + " " + array[0]

        elif kind == 3:
            # FOR DIHEDRALS, IF COLUMN 3 HAS LOWER PRIORITY THAN COLUMN 2, THE ORDER FOR THE ENTIRE DIHEDRAL IS REVERSED
            # FOR DIHEDRALS, IF COLUMNS 2 & 3 HAVE THE SAME PRIORITY, COLUMS 1 & 4 ARE CONSIDERED INSTEAD.
            #   IF 4 HAS LOWER PRIORITY THAN 1, THE ORDER FOR THE ENTIRE DIHEDRAL IS REVERSED
            array = atoms_combi.split(" ")
            atom_priority_i = int(hAtomPriority[array[0]])
            atom_priority_j = int(hAtomPriority[array[1]])
            atom_priority_k = int(hAtomPriority[array[2]])
            atom_priority_l = int(hAtomPriority[array[3]])
            if atom_priority_j < atom_priority_k:
                return atoms_combi
            else:
                if atom_priority_k < atom_priority_j:
                    return array[3] + " " + array[2] + " " + array[1] + " " + array[0]
                elif atom_priority_k == atom_priority_j:
                    if atom_priority_i <= atom_priority_l:
                        return atoms_combi
                    else:
                        return array[3] + " " + array[2] + " " + array[1] + " " + array[0]

        elif kind == 4:
            # - COLUMN 1 IS ALWAYS THE CENTRAL ATOM
            # - IF 2 OF THE SUBSTITUENTS HAVE IDENTICAL TYPES, THESE SHOULD
            #   BE IN COLUMNS 2 & 3 (BUT THEY CANNOT BE MOVED AROUND
            #   WITHOUT RE-OPTIMIZING THE PARAMETER)
            # - IF THE SUBSTITUENTS ARE ALL DIFFERENT, COLUMNS 2, 3 & 4
            #   SHOULD BE SORTED BY INCREASING PRIORITY. COLUMNS 2 AND 3
            #   CAN BE SWAPPED WITHOUT CHANGING THE PARAMETER BUT OTHER
            #   PERMUTATIONS MANDATE RE-OPTIMIZATION
            array = atoms_combi.split(" ")
            atom_priority_i = int(hAtomPriority[array[0]])
            atom_priority_j = int(hAtomPriority[array[1]])
            atom_priority_k = int(hAtomPriority[array[2]])
            atom_priority_l = int(hAtomPriority[array[3]])

            if atom_priority_j == atom_priority_k:
                return array[0] + " " + array[1] + " " + array[2] + " " + array[3]
            elif atom_priority_j == atom_priority_l:
                return array[0] + " " + array[1] + " " + array[3] + " " + array[2]
            elif atom_priority_k == atom_priority_l:
                return array[0] + " " + array[2] + " " + array[3] + " " + array[1]
            else:
                index = []
                index.append(1)
                index.append(2)
                index.append(3)

                priority = []
                priority.append(atom_priority_j)
                priority.append(atom_priority_k)
                priority.append(atom_priority_l)

                for i in range(3):
                    for j in range(i+1, 3):
                        if priority[i] > priority[j]:
                            temp_index = index[j]
                            index[j] = index[i]
                            index[i] = temp_index

                            temp_priority = priority[j]
                            priority[j] = priority[i]
                            priority[i] = temp_priority

                return array[0] + " " + array[index[0]] + " " + array[index[1]] + " " + array[index[2]]

        return atoms_combi



    def orderAtomsBonds(self, hAtomPriority, mol):
        atom_priority_i = 0
        atom_priority_j = 0

        for i in range(len(mol.bonds)):
            atom_priority_i = int(hAtomPriority[mol.atoms[mol.bonds[i].i].atomType])
            atom_priority_j = int(hAtomPriority[mol.atoms[mol.bonds[i].j].atomType])

            if atom_priority_i > atom_priority_j:
                temp = mol.bonds[i].i
                mol.bonds[i].i = mol.bonds[i].j
                mol.bonds[i].j = temp

    def orderAtomsAngles(self, hAtomPriority, mol):
        atom_priority_i = 0
        atom_priority_j = 0
        atom_priority_k = 0

        for i in range(len(mol.angles)):
            array = mol.angles[i].linked_atoms.split( "-")
            atom_priority_i = int(hAtomPriority[mol.atoms[int(array[0])].atomType])
            atom_priority_k = int(hAtomPriority[mol.atoms[int(array[2])].atomType])

            if atom_priority_i > atom_priority_k:
                mol.angles[i].linked_atoms = array[2] + "-" + array[1] + "-" + array[0]

    def setFormalCharge(self, mol):
        num_atoms = len(mol.atoms)
        atom = Atom()
        i = 0
        while i < num_atoms:
            atom = mol.atoms[i]
            #  -1
            #  heavy atoms linked to OG2D2, OG312, SG302, OG2P1
            #  PG1 linked to OG2P1
            if atom.atomType.lower() == "OG2D2".lower() or atom.atomType.lower() == "OG312".lower() or atom.atomType.lower() == "SG302".lower() or atom.atomType.lower() == "OG2P1".lower():
                while j < 5:
                    if atom.linkage[j] != -1:
                        if mol.atoms[atom.linkage[j]].element[0] == 'C' or mol.atoms[atom.linkage[j]].element[0] == 'S':
                            mol.atoms[atom.linkage[j]].formal_charge = -1
                    j += 1
            if atom.atomType.lower() == "OG2P1".lower():
                if mol.atoms[atom.linkage[0]].atomType.lower() == "PG1".lower():
                    #  -1
                    #  PG1 linked to OG2P1
                    mol.atoms[atom.linkage[0]].formal_charge = -1
                elif mol.atoms[atom.linkage[0]].atomType.lower() == "PG2".lower():
                    #  -2
                    #  PG2 linked to to OG2P1
                    mol.atoms[atom.linkage[0]].formal_charge = -2
            #  1
            #  heavy atoms linked to NG2P1, NG3P3, NG2R52, NG3P2, NG3P0, NG2R61, NG3P1
            if atom.atomType.lower() == "NG2P1".lower() or atom.atomType.lower() == "NG3P3".lower() or atom.atomType.lower() == "NG2R52".lower() or atom.atomType.lower() == "NG3P2".lower() or atom.atomType.lower() == "NG3P0".lower() or atom.atomType.lower() == "NG2R61".lower() or atom.atomType.lower() == "NG3P1".lower():
                while j < 5:
                    if atom.linkage[j] != -1:
                        if mol.atoms[atom.linkage[j]].element[0] == 'C':
                            mol.atoms[atom.linkage[j]].formal_charge = 1
                    j += 1
            i += 1

    def removeZeroIncrementBond(self, mol):
        #  The charge increment for a bond between two identical atom types is zero by definition.
        i = 0
        while i < len(mol.bonds):
            if mol.atoms[mol.bonds[i].i].atomType.lower() == mol.atoms[mol.bonds[i.lower().j].atomType]:
                del mol.bonds[i]
                i -= 1
            elif mol.atoms[mol.bonds[i].i].atomType.lower() == "CG2D1O".lower() and mol.atoms[mol.bonds[i].j].atomType.lower() == "CG2D2O".lower():
                del mol.bonds[i]
                i -= 1
            elif mol.atoms[mol.bonds[i].i].atomType.lower() == "CG2D2O".lower() and mol.atoms[mol.bonds[i].j].atomType.lower() == "CG2D1O".lower():
                del mol.bonds[i]
                i -= 1
            elif mol.atoms[mol.bonds[i].i].atomType.lower() == "CG2DC1".lower() and mol.atoms[mol.bonds[i].j].atomType.lower() == "CG2DC2".lower():
                del mol.bonds[i]
                i -= 1
            elif mol.atoms[mol.bonds[i].i].atomType.lower() == "CG2DC2".lower() and mol.atoms[mol.bonds[i].j].atomType.lower() == "CG2DC1".lower():
                del mol.bonds[i]
                i -= 1
            elif mol.atoms[mol.bonds[i].i].atomType.lower() == "CG25C1".lower() and mol.atoms[mol.bonds[i].j].atomType.lower() == "CG25C2".lower():
                del mol.bonds[i]
                i -= 1
            elif mol.atoms[mol.bonds[i].i].atomType.lower() == "CG25C2".lower() and mol.atoms[mol.bonds[i].j].atomType.lower() == "CG25C1".lower():
                del mol.bonds[i]
                i -= 1
            elif mol.atoms[mol.bonds[i].i].atomType.lower() == "CG251O".lower() and mol.atoms[mol.bonds[i].j].atomType.lower() == "CG252O".lower():
                del mol.bonds[i]
                i -= 1
            elif mol.atoms[mol.bonds[i].i].atomType.lower() == "CG252O".lower() and mol.atoms[mol.bonds[i].j].atomType.lower() == "CG251O".lower():
                del mol.bonds[i]
                i -= 1
            i += 1

    def setChargeIncrementBond(self, mol, hChargeIncrement):
        value = None
        atom_combi = ""
        i = 0
        while i < len(mol.bonds):
            atom_combi = mol.atoms[mol.bonds[i].i].atomType + " " + mol.atoms[mol.bonds[i].j].atomType
            value = hChargeIncrement.get(atom_combi)
            if value != None:
                mol.bonds[i].increment = Double.parseDouble(value)
            else:
                print "No available bond charge increment: " + atom_combi
            i += 1

    def assignAtomChargesUsingBondIncrement(self, mol):
        a = 0
        while a < len(mol.atoms):
            mol.atoms[a].charge_calc = mol.atoms[a].formal_charge
            a += 1
        b = 0
        while b < len(mol.bonds):
            mol.atoms[mol.bonds[b].i].charge_calc -= mol.bonds[b].increment
            mol.atoms[mol.bonds[b].j].charge_calc += mol.bonds[b].increment
            b += 1

    #  =============================================================================
    #  Algorithm is based on Hanser et al. J. Chem. Inf. Comput. Sci. 1996, 36, 1146
    #  A new algorithm for exhaustive ring perception in a molecular graph
    #  =============================================================================
    def detectRings(self, mol):
        atom = Atom()
        edge = Edge()
        frag = Fragment()
        vEdges = []
        vRings = []
        #  *** generate P-Graph ***
        vEdges = self.getInitialEdgesUsingBonds(mol)
        while self.getNumRemainingVertex(mol.atoms) != 0:
            vertex = self.chooseVertex( mol.atoms )
            self.clearEdges(vEdges, 7)
            self.remove(vertex, mol.atoms, vEdges, vRings)
        self.removeNonRings(mol.atoms, vRings)
        self.removeRedundantRings(vRings)
        self.keepSmallestRings(vRings)
        return vRings

    def setRingAtoms(self, mol, vRings):
        array = []
        index = int()
        i = 0
        while i < len(vRings):
            array = vRings[i].path.split(" ")
            while len(array):
                index = int(array[j])
                mol.atoms[index].isRingAtom = True
                while k < 3:
                    if mol.atoms[index].ringIndex[k] == -1:
                        mol.atoms[index].ringIndex[k] = i
                        len(array)
                        break
                    k += 1
                j += 1
            i += 1

    def setRingExtraProperties(self, mol, vRings, aRingAromacity, aRingHasCabonyl):
        array = []
        index = int()
        is_pyrrole = False
        #  aromatic
        is_neutral_pyridine = False
        #  aromatic
        is_protonated_pyridine = False
        is_protonated_pyrimidine = False
        #  aromatic
        is_furan = False
        #  aromatic
        is_furanose = False
        i = 0
        while i < len(vRings):
            is_pyrrole = False
            is_neutral_pyridine = False
            is_protonated_pyridine = False
            is_protonated_pyrimidine = False
            array = vRings[i].path.split(" ")
            # print  vRings[i).path ;
            #  *** check planarity ***
            if len(array):
                if isPlanarRing(mol, array):
                    aRingAromacity[i] = True
                    while len(array):
                        index = int(array[j])
                        mol.atoms[index].isAromatic = True
                        if mol.atoms[index].isCarbonyl:
                            aRingHasCabonyl[i] = True
                        j += 1
            if isPyrrole(mol, array):
                is_pyrrole = True
            if isNeutralPyridine(mol, array):
                is_neutral_pyridine = True
            if isProtonatedPyridine(mol, array):
                is_protonated_pyridine = True
            if isProtonatedPyrimidine(mol, array):
                is_protonated_pyrimidine = True
            if isFuran(mol, array):
                is_furan = True
            if isFuranose(mol, array):
                is_furanose = True
            while len(array):
                index = int(array[j])
                mol.atoms[index].isPyrrole = is_pyrrole
                mol.atoms[index].isNeutralPyridine = is_neutral_pyridine
                mol.atoms[index].isProtonatedPyridine = is_protonated_pyridine
                mol.atoms[index].isProtonatedPyrimidine = is_protonated_pyrimidine
                mol.atoms[index].isFuran = is_furan
                mol.atoms[index].isFuranose = is_furanose
                j += 1
            i += 1

    def isPlanarRing(self, mol, array):
        atom = Atom()
        i = 0
        while len(array):
            atom = mol.atoms[int(array[i])]
            if atom.num_linkages == 4:
                return False
            i += 1
        return True

    def isPyrrole(self, mol, array):
        atom = Atom()
        num_carbon_atoms = 0
        num_nitrogen_atoms = 0
        if len(array):
            while i < 5:
                atom = mol.atoms[int(array[i])]
                if not atom.isAromatic:
                    return False
                if atom.element.lower() == "C".lower():
                    num_carbon_atoms += 1
                if atom.element.lower() == "N".lower():
                    num_nitrogen_atoms += 1
                i += 1
            if num_carbon_atoms == 4 and num_nitrogen_atoms == 1:
                return True
        return False

    def isNeutralPyridine(self, mol, array):
        atom = Atom()
        num_carbon_atoms = 0
        num_nitrogen_atoms = 0
        num_protonated_N = 0
        if len(array):
            while i < 6:
                atom = mol.atoms[int(array[i])]
                if not atom.isAromatic:
                    return False
                if atom.element.lower() == "C".lower():
                    num_carbon_atoms += 1
                if atom.element.lower() == "N".lower():
                    num_nitrogen_atoms += 1
                    if atom.isProtonatedNitrogen:
                        num_protonated_N += 1
                i += 1
            if num_carbon_atoms == 5 and num_nitrogen_atoms == 1 and num_protonated_N == 0:
                return True
        return False

    def isProtonatedPyridine(self, mol, array):
        atom = Atom()
        num_carbon_atoms = 0
        num_nitrogen_atoms = 0
        num_protonated_N = 0
        if len(array):
            while i < 6:
                atom = mol.atoms[int(array[i])]
                if not atom.isAromatic:
                    return False
                if atom.element.lower() == "C".lower():
                    num_carbon_atoms += 1
                if atom.element.lower() == "N".lower():
                    num_nitrogen_atoms += 1
                    if atom.isProtonatedNitrogen:
                        num_protonated_N += 1
                i += 1
            if num_carbon_atoms == 5 and num_nitrogen_atoms == 1 and num_protonated_N != 0:
                return True
        return False

    def isProtonatedPyrimidine(self, mol, array):
        atom = Atom()
        num_carbon_atoms = 0
        num_nitrogen_atoms = 0
        num_protonated_N = 0
        if len(array):
            while i < 6:
                atom = mol.atoms[int(array[i])]
                if not atom.isAromatic:
                    return False
                if atom.element.lower() == "C".lower():
                    num_carbon_atoms += 1
                if atom.element.lower() == "N".lower():
                    num_nitrogen_atoms += 1
                    if atom.isProtonatedNitrogen:
                        num_protonated_N += 1
                i += 1
            if num_carbon_atoms == 4 and num_nitrogen_atoms == 2 and num_protonated_N != 0:
                return True
        return False

    def isFuran(self, mol, array):
        atom = Atom()
        num_carbon_atoms = 0
        num_nitrogen_atoms = 0
        num_oxygen_atoms = 0
        isAromatic = True
        if len(array):
            while i < 5:
                atom = mol.atoms[int(array[i])]
                if not atom.isAromatic:
                    isAromatic = False
                if atom.element.lower() == "C".lower():
                    num_carbon_atoms += 1
                if atom.element.lower() == "N".lower():
                    num_nitrogen_atoms += 1
                if atom.element.lower() == "O".lower():
                    num_oxygen_atoms += 1
                i += 1
            if num_carbon_atoms == 4 and num_oxygen_atoms == 1 and isAromatic:
                return True
            elif num_carbon_atoms == 3 and num_nitrogen_atoms == 1 and num_oxygen_atoms == 1 and isAromatic:
                return True
            elif num_carbon_atoms == 2 and num_nitrogen_atoms == 2 and num_oxygen_atoms == 1 and isAromatic:
                return True
        return False

    def isFuranose(self, mol, array):
        atom = Atom()
        index = int()
        num_carbon_atoms = 0
        num_nitrogen_atoms = 0
        num_oxygen_atoms = 0
        isAromatic = True
        if len(array):
            while i < 5:
                index = int(array[i])
                atom = mol.atoms[index]
                if not atom.isAromatic:
                    isAromatic = False
                if atom.element.lower() == "C".lower():
                    num_carbon_atoms += 1
                if atom.element.lower() == "N".lower():
                    num_nitrogen_atoms += 1
                if atom.element.lower() == "O".lower():
                    num_oxygen_atoms += 1
                i += 1
            if num_carbon_atoms == 4 and num_oxygen_atoms == 1 and not isAromatic:
                return True
            elif num_carbon_atoms == 3 and num_nitrogen_atoms == 1 and num_oxygen_atoms == 1 and not isAromatic:
                return True
            elif num_carbon_atoms == 3 and num_oxygen_atoms == 2:
                return True
        return False

    def setNumHydrogenAtoms(self, mol):
        num_hydrogen_atoms = 0
        i = 0
        while i < len(mol.atoms):
            num_hydrogen_atoms = 0
            j = 0
            while j < 5:
                if mol.atoms[i].linkage[j] != -1:
                    if mol.atoms[mol.atoms[i].linkage[j]].element.lower() == "H".lower():
                        num_hydrogen_atoms += 1
                j += 1
            mol.atoms[i].numHydrogenAtoms = num_hydrogen_atoms
            i += 1

    def setNumCarbonAtoms(self, mol):
        num_carbon_atoms = 0
        i = 0
        while i < len(mol.atoms):
            num_carbon_atoms = 0
            j = 0
            while j < 5:
                if mol.atoms[i].linkage[j] != -1:
                    if mol.atoms[mol.atoms[i].linkage[j]].element.lower() == "C".lower():
                        num_carbon_atoms += 1
                j += 1
            mol.atoms[i].numCarbonAtoms = num_carbon_atoms
            i += 1

    def setNumNitrogenAtoms(self, mol):
        num_nigrogen_atoms = 0
        i = 0
        while i < len(mol.atoms):
            num_nigrogen_atoms = 0
            j = 0
            while j < 5:
                if mol.atoms[i].linkage[j] != -1:
                    if mol.atoms[mol.atoms[i].linkage[j]].element.lower() == "N".lower():
                        num_nigrogen_atoms += 1
                j += 1
            mol.atoms[i].numNitrogenAtoms = num_nigrogen_atoms
            i += 1

    def setNumOxygenAtoms(self, mol):
        num_oxygen_atoms = 0
        i = 0
        while i < len(mol.atoms):
            num_oxygen_atoms = 0
            j = 0
            while j < 5:
                if mol.atoms[i].linkage[j] != -1:
                    if mol.atoms[mol.atoms[i].linkage[j]].element.lower() == "O".lower():
                        num_oxygen_atoms += 1
                j += 1
            mol.atoms[i].numOxygenAtoms = num_oxygen_atoms
            i += 1

    def setNumSulfurAtoms(self, mol):
        num_sulfur_atoms = 0
        i = 0
        while i < len(mol.atoms):
            num_sulfur_atoms = 0
            j = 0
            while j < 5:
                if mol.atoms[i].linkage[j] != -1:
                    if mol.atoms[mol.atoms[i].linkage[j]].element.lower() == "S".lower():
                        num_sulfur_atoms += 1
                j += 1
            mol.atoms[i].numSulfurAtoms = num_sulfur_atoms
            i += 1

    def setMethylAtoms(self, mol):
        i = 0
        while i < len(mol.atoms):
            if mol.atoms[i].element.lower() == "C".lower() and mol.atoms[i].numHydrogenAtoms == 3:
                mol.atoms[i].isMethyl = True
            i += 1

    def isAtomInFiveMemberedRing(self, mol, atm_index):
        atom = mol.atoms[atm_index]
        i = 0
        while i < 3:
            if atom.numRingAtoms[i] == 5:
                return True
            i += 1
        return False

    def isAtomInSixMemberedRing(self, mol, atm_index):
        atom = mol.atoms[atm_index]
        i = 0
        while i < 3:
            if atom.numRingAtoms[i] == 6:
                return True
            i += 1
        return False

    #  ==============================================
    #              Assign Bond Order
    #  ==============================================
    #  element / sum of bond orders
    #  H: 1
    #  C: 4
    #  N: 3 (0), 4 (+)
    #  O: 2 (0), 1 (-1)
    #  S: 2 (0), 4 (0), 6 (0), 1 (-1), 4 (-1), 7 (-1), 6 (-1)
    #  CL, BR, I, F: 1
    #  P: 5 (0), 6 (-1), 7 (-2)
    def setBondOrders(self, mol):
        atom = Atom()
        atom_linked = Atom()
        #  assign bond orders using bond information in mol2
        i = 0
        while i < len(mol.bonds):
            j = 0
            while j < 5:
                if mol.atoms[mol.bonds[i].i].linkage[j] == mol.bonds[i].j:
                    if mol.bonds[i].bondType[0] == '1':
                        mol.atoms[mol.bonds[i].i].bondOrder[j] = 1
                    elif mol.bonds[i].bondType[0] == '2':
                        mol.atoms[mol.bonds[i].i].bondOrder[j] = 2
                    elif mol.bonds[i].bondType[0] == '3':
                        mol.atoms[mol.bonds[i].i].bondOrder[j] = 3
                    elif mol.bonds[i].bondType.lower() == "am".lower():
                        mol.atoms[mol.bonds[i].i].bondOrder[j] = 1
                    break
                j += 1
            while j < 5:
                if mol.atoms[mol.bonds[i].j].linkage[j] == mol.bonds[i].i:
                    if mol.bonds[i].bondType[0] == '1':
                        mol.atoms[mol.bonds[i].j].bondOrder[j] = 1
                    elif mol.bonds[i].bondType[0] == '2':
                        mol.atoms[mol.bonds[i].j].bondOrder[j] = 2
                    elif mol.bonds[i].bondType[0] == '3':
                        mol.atoms[mol.bonds[i].j].bondOrder[j] = 3
                    elif mol.bonds[i].bondType.lower() == "am".lower():
                        mol.atoms[mol.bonds[i].j].bondOrder[j] = 1
                    break
                j += 1
            i += 1
        i = 0
        while i < len(mol.atoms):
            atom = mol.atoms[i]
            if atom.isAromatic and atom.element.lower() == "C".lower() and self.hasDoubleBondedOxygen(mol, i) and atom.num_linkages == 3:
                while j < 3:
                    if not mol.atoms[atom.linkage[j]].element.lower() == "O".lower():
                        atom.bondOrder[j] = 1
                        while k < 5:
                            atom_linked = mol.atoms[atom.linkage[j]]
                            if atom_linked.linkage[k] == i:
                                atom_linked.bondOrder[k] = 1
                                break
                            k += 1
                    j += 1
            i += 1
        i = 0
        while i < len(mol.atoms):
            atom = mol.atoms[i]
            if atom.isAromatic and atom.element.lower() == "C".lower() and atom.num_linkages == 3 and getNumUnassignedBondOrders(atom) == 1:
                while j < 3:
                    if atom.bondOrder[j] == 0:
                        atom.bondOrder[j] = 2
                        while k < 5:
                            atom_linked = mol.atoms[atom.linkage[j]]
                            if atom_linked.linkage[k] == i:
                                atom_linked.bondOrder[k] = 2
                                break
                            k += 1
                        break
                    j += 1
            i += 1
        i = 0
        while i < len(mol.atoms):
            atom = mol.atoms[i]
            if atom.isAromatic and atom.element.lower() == "C".lower() and atom.num_linkages == 3 and getNumUnassignedBondOrders(atom) == 2:
                #  double bond for one bond
                while j < 3:
                    if atom.bondOrder[j] == 0:
                        atom.bondOrder[j] = 2
                        while k < 5:
                            atom_linked = mol.atoms[atom.linkage[j]]
                            if atom_linked.linkage[k] == i:
                                atom_linked.bondOrder[k] = 2
                                break
                            k += 1
                        break
                    j += 1
                #  single bond for the other bond
                while j < 3:
                    if atom.bondOrder[j] == 0:
                        atom.bondOrder[j] = 1
                        while k < 5:
                            atom_linked = mol.atoms[atom.linkage[j]]
                            if atom_linked.linkage[k] == i:
                                atom_linked.bondOrder[k] = 1
                                break
                            k += 1
                        break
                    j += 1
            i += 1

    #  ==============================================
    def getNumUnassignedBondOrders(self, atom):
        num_unassigned_bond_orders = 0
        i = 0
        while i < atom.num_linkages:
            if atom.bondOrder[i] == 0:
                num_unassigned_bond_orders += 1
            i += 1
        return num_unassigned_bond_orders

    def setTotalBondOrder(self, mol):
        atom = Atom()
        num_atoms = len(mol.atoms)
        total_bond_order = 0
        i = 0
        j = 0
        while i < num_atoms:
            atom = mol.atoms[i]
            total_bond_order = 0
            while j < 5:
                total_bond_order += atom.bondOrder[j]
                j += 1
            mol.atoms[i].totalBondOrder = total_bond_order
            i += 1

    def setAngleList(self, mol):
        atmgrp = LinkedAtomGroup()
        array = []
        i = 0
        while i < len(mol.bonds):
            for j in range(i + 1, len(mol.bonds)):
                if mol.bonds[i].i == mol.bonds[j].j:
                    atmgrp = LinkedAtomGroup()
                    atmgrp.i = mol.bonds[j].i
                    atmgrp.j = mol.bonds[i].j
                    atmgrp.linked_atoms = "" + str(mol.bonds[j].i) + "-" + str(mol.bonds[i].i) + "-" + str(mol.bonds[i].j)
                    atmgrp.linked_bonds = "" + str(j) + "-" + str(i)
                    mol.angles.append(atmgrp)
                elif mol.bonds[i].j == mol.bonds[j].i:
                    atmgrp = LinkedAtomGroup()
                    atmgrp.i = mol.bonds[i].i
                    atmgrp.j = mol.bonds[j].j
                    atmgrp.linked_atoms = "" + str(mol.bonds[i].i) + "-" + str(mol.bonds[i].j) + "-" + str(mol.bonds[j].j)
                    atmgrp.linked_bonds = "" + str(i) + "-" + str(j)
                    mol.angles.append(atmgrp)
                elif mol.bonds[i].i == mol.bonds[j].i:
                    atmgrp = LinkedAtomGroup()
                    atmgrp.i = mol.bonds[j].j
                    atmgrp.j = mol.bonds[i].j
                    atmgrp.linked_atoms = "" + str(mol.bonds[j].j) + "-" + str(mol.bonds[i].i) + "-" + str(mol.bonds[i].j)
                    atmgrp.linked_bonds = "" + str(j) + "-" + str(i)
                    mol.angles.append(atmgrp)
                elif mol.bonds[i].j == mol.bonds[j].j:
                    atmgrp = LinkedAtomGroup()
                    atmgrp.i = mol.bonds[i].i
                    atmgrp.j = mol.bonds[j].i
                    atmgrp.linked_atoms = "" + str(mol.bonds[i].i) + "-" + str(mol.bonds[i].j) + "-" + str(mol.bonds[j].i)
                    atmgrp.linked_bonds = "" + str(i) + "-" + str(j)
                    mol.angles.append(atmgrp)
            i += 1
        #
        #       print  len(mol.angles) ;
        #       for( int i = 0 ; i < len(mol.angles) ; i++ ){
        #           for( int j = i+1 ; j < len(mol.angles) ; j++ ){
        #               array = mol.angles[j).linked_atoms.split( "-" );
        #               if( mol.angles[i).linked_atoms.lower() ==  mol.angles[j.lower().linked_atoms ) ){
        #                   mol.angles.removeElementAt(j);
        #                   j--;
        #               }
        #           }
        #       }
        #       print  len(vAngle) ;
        #

    def setDihderalList(self, mol):
        atmgrp = LinkedAtomGroup()
        array = []
        i = 0
        while i < len(mol.angles):
            array = mol.angles[i].linked_bonds.split("-")
            j = 0
            while j < len(mol.bonds):
                if j != int(array[0]) and j != int(array[1]):
                    if mol.angles[i].i == mol.bonds[j].j:
                        atmgrp = LinkedAtomGroup()
                        atmgrp.i = mol.bonds[j].i
                        atmgrp.j = mol.angles[i].j
                        atmgrp.linked_atoms = str(mol.bonds[j].i) + "-" + str(mol.angles[i].linked_atoms)
                        atmgrp.linked_bonds = str(j) + "-" + str(mol.angles[i].linked_bonds)
                        mol.dihedrals.append(atmgrp)
                    elif mol.angles[i].j == mol.bonds[j].i:
                        atmgrp = LinkedAtomGroup()
                        atmgrp.i = mol.angles[i].i
                        atmgrp.j = mol.bonds[j].j
                        atmgrp.linked_atoms = str(mol.angles[i].linked_atoms) + "-" + str(mol.bonds[j].j)
                        atmgrp.linked_bonds = str(mol.angles[i].linked_bonds) + "-" + str(j)
                        mol.dihedrals.append(atmgrp)
                    elif mol.angles[i].i == mol.bonds[j].i:
                        atmgrp = LinkedAtomGroup()
                        atmgrp.i = mol.bonds[j].j
                        atmgrp.j = mol.angles[i].j
                        atmgrp.linked_atoms = str(mol.bonds[j].j) + "-" + str(mol.angles[i].linked_atoms)
                        atmgrp.linked_bonds = str(j) + "-" + str(mol.angles[i].linked_bonds)
                        mol.dihedrals.append(atmgrp)
                    elif mol.angles[i].j == mol.bonds[j].j:
                        atmgrp = LinkedAtomGroup()
                        atmgrp.i = mol.angles[i].i
                        atmgrp.j = mol.bonds[j].i
                        atmgrp.linked_atoms = str(mol.angles[i].linked_atoms) + "-" + str(mol.bonds[j].i)
                        atmgrp.linked_bonds = str(mol.angles[i].linked_bonds) + "-" + str(j)
                        mol.dihedrals.append(atmgrp)
                j += 1
            i += 1
        # print  len(vDihedral) ;
        i = 0
        while i < len(mol.dihedrals):
            while j < len(mol.dihedrals):
                if mol.dihedrals[i].linked_atoms.lower() == mol.dihedrals[j].linked_atoms.lower():
                    del mol.dihedrals[j]
                    j -= 1
                j += 1
            i += 1

    def setImproperList(self, mol, aRingAromacity):
        num_atoms = len(mol.atoms)
        atom = Atom()
        atom_type = str()
        atmgrp = LinkedAtomGroup()
        add_to_list = False
        i = 0
        while i < num_atoms:
            add_to_list = False
            atom = mol.atoms[i]
            atom_type = atom.atomType
            if atom.element.lower() == "C".lower() and atom.num_linkages == 3:
                if not atom.isRingAtom and self.getNumberOfLinkedRingAtoms(mol, i) <= 1:
                    if atom.isCarbonyl and atom.numCarbonAtoms == 1 and atom.numOxygenAtoms == 2:
                        add_to_list = True
                    elif atom.isCarbonyl and atom.numCarbonAtoms == 1 and atom.numOxygenAtoms == 2:
                        add_to_list = True
                    elif atom.isCarbonyl and atom.numNitrogenAtoms == 1 and atom.numOxygenAtoms == 1 and atom.numHydrogenAtoms == 1:
                        add_to_list = True
                    elif atom.isCarbonyl and atom.numOxygenAtoms == 3:
                        add_to_list = True
                    elif atom.isCarbonyl and atom.numNitrogenAtoms == 1 and atom.numOxygenAtoms == 2:
                        add_to_list = True
                    elif atom.isImineCarbon and atom.numNitrogenAtoms == 3:
                        add_to_list = True
                    elif atom.isCarbonyl and atom.numCarbonAtoms == 1 and atom.numNitrogenAtoms == 1 and atom.numOxygenAtoms == 1:
                        add_to_list = True
                    elif atom.isCarbonyl and atom.numCarbonAtoms == 1 and atom.numOxygenAtoms == 1 and atom.numHydrogenAtoms == 1:
                        add_to_list = True
                    elif atom.isCarbonyl and atom.numCarbonAtoms == 2 and atom.numOxygenAtoms == 1:
                        add_to_list = True
                    elif atom.isCarbonyl and atom.numOxygenAtoms == 2 and atom.numHydrogenAtoms == 1:
                        add_to_list = True
                    elif atom.isCarbonyl and atom.numSulfurAtoms == 3:
                        add_to_list = True
                    elif atom.isImineCarbon and atom.numCarbonAtoms == 1 and atom.numNitrogenAtoms == 1 and atom.numHydrogenAtoms == 1:
                        add_to_list = True
                    elif atom.isImineCarbon and atom.numCarbonAtoms == 1 and atom.numNitrogenAtoms == 2:
                        add_to_list = True
                    elif hasDoubleBondedCarbon(mol, i) and atom.numNitrogenAtoms == 1 and atom.numCarbonAtoms == 1 and atom.numHydrogenAtoms == 1:
                        add_to_list = True
                    elif hasDoubleBondedCarbon(mol, i) and atom.numOxygenAtoms == 1 and atom.numCarbonAtoms == 1 and atom.numHydrogenAtoms == 1:
                        add_to_list = True
                    elif hasDoubleBondedCarbon(mol, i) and atom.numCarbonAtoms == 2 and atom.numNitrogenAtoms == 1:
                        add_to_list = True
                elif atom.isRingAtom and self.getNumberOfLinkedRingAtoms(mol, i) == 2:
                    if atom.isImineCarbon and getFirstExocylicAtom(mol, i).lower() == "N".lower() and atom.numCarbonAtoms == 1 and atom.numNitrogenAtoms == 2:
                        add_to_list = True
                    elif atom.isImineCarbon and atom.numNitrogenAtoms == 3:
                        add_to_list = True
                    elif atom.isCarbonyl and atom.numNitrogenAtoms == 1 and atom.numOxygenAtoms == 1 and atom.numCarbonAtoms == 1:
                        add_to_list = True
                    elif atom.isCarbonyl and atom.numNitrogenAtoms == 2 and atom.numOxygenAtoms == 1:
                        add_to_list = True
                    elif atom.isCarbonyl and atom.numCarbonAtoms == 2 and atom.numOxygenAtoms == 1:
                        add_to_list = True
                    elif hasDoubleBondedCarbon(mol, i) and not atom.isAromatic and atom.numNitrogenAtoms == 1 and atom.numCarbonAtoms == 1 and atom.numHydrogenAtoms == 1:
                        add_to_list = True
                    elif atom.isCarbonyl and atom.numNitrogenAtoms == 1 and atom.numSulfurAtoms == 2:
                        add_to_list = True
                    elif atom.isCarbonyl and atom.numCarbonAtoms == 1 and atom.numSulfurAtoms == 1 and atom.numNitrogenAtoms == 1:
                        add_to_list = True
                    elif atom.isCarbonyl and atom.numNitrogenAtoms == 2 and atom.numSulfurAtoms == 1:
                        add_to_list = True
                    elif atom.isImineCarbon and atom.numNitrogenAtoms == 2 and atom.numSulfurAtoms == 1:
                        add_to_list = True
                    elif atom.isCarbonyl and atom.numCarbonAtoms == 1 and atom.numOxygenAtoms == 2:
                        add_to_list = True
                    elif atom.isCarbonyl and atom.numNitrogenAtoms == 1 and atom.numOxygenAtoms == 2:
                        add_to_list = True
                    elif atom.isImineCarbon and atom.numNitrogenAtoms == 2 and atom.numOxygenAtoms == 1:
                        add_to_list = True
            elif atom.element.lower() == "N".lower() and atom.num_linkages == 3:
                if not atom.isRingAtom and self.getNumberOfLinkedRingAtoms(mol, i) <= 1:
                    if self.hasDoubleBondedOxygen(mol, i) and atom.numCarbonAtoms == 1 and atom.numOxygenAtoms == 2:
                        add_to_list = True
                    elif self.getNumberOfLinkedRingAtoms(mol, i) == 1 and atom.numCarbonAtoms == 1 and atom.numHydrogenAtoms == 2:
                        add_to_list = True
            if add_to_list:
                atmgrp = LinkedAtomGroup()
                atmgrp.linked_atoms = "" + str(i)
                j = 0
                while j < 3:
                    atmgrp.linked_atoms += "-" + str(atom.linkage[j])
                    j += 1
                mol.impropers.append(atmgrp)
            i += 1


    def getNumberOfLinkedRingAtoms(self, mol, atm_index):
        atom = mol.atoms[atm_index]
        num_ring_atoms = 0
        i = 0
        while i < atom.num_linkages:
            if mol.atoms[atom.linkage[i]].isRingAtom:
                num_ring_atoms += 1
            i += 1
        return num_ring_atoms

    def getFirstExocylicAtom(self, mol, atm_index):
        atom = mol.atoms[atm_index]
        i = 0
        while i < atom.num_linkages:
            if not mol.atoms[atom.linkage[i]].isRingAtom:
                return mol.atoms[atom.linkage[i]].element
            i += 1
        return "X"

    def sortAtoms(self, atm_grp, kind):
        #  ORDER OF PREFERENCE FOR SORTING PARAMETERS:
        #  C < N < O < P < S < HALOGENS (LOW TO HIGH Z) < MISC. (BY Z) < H
        #  ATOMS TYPES WITHIN THE SAME ELEMENT ARE SORTED ALPHABETICALLY
        atoms = ""
        if kind == 1:
            """#  BOND"""
            #  THE LOWEST PRIORITY ATOM ALWAYS COMES FIRST
        elif kind == 2:
            """#  ANGLE"""
        elif kind == 3:
            """#  DIHEDRAL"""
        elif kind == 4:
            """#  IMPROPER"""
        return atoms

    def setConjugatedAtoms(self, mol):
        # C=C-C=C: 212
        # C=C-C=O: 212
        # N=C-C=C: 212
        # C=C(al)-C=C(ar): 21ar

        # 1 = single
        # 2 = double
        # ar = aromatic
        consecutive_elements = ""
        consecutive_bonds = ""

        for i in range(len(mol.dihedrals)):
            array = mol.dihedrals[i].linked_atoms.split("-")
            consecutive_elements = ""
            consecutive_bonds = ""

            for j in range(4):
                consecutive_elements = consecutive_elements + mol.atoms[int(array[j])].element

            consecutive_elements = consecutive_elements.strip()

            array = mol.dihedrals[i].linked_bonds.split("-")

            for j in range(3):
                consecutive_bonds = consecutive_bonds + mol.bonds[int(array[j])].bondType

            if (consecutive_elements.lower() == "CCCC".lower()) or (consecutive_elements.lower() == "CCC0".lower()) or (consecutive_elements.lower() == "0CCC".lower())\
                or (consecutive_elements.lower() == "NCCC".lower()) or (consecutive_elements.lower() == "CCCN".lower()) :
                if (consecutive_bonds.lower() == "212") or (consecutive_bonds.lower() == "21ar") or (consecutive_bonds.lower() == "ar12"):
                    array = mol.dihedrals[i].linked_atoms.split("-")

                    for j in range(4):
                        mol.atoms[int(array[j])].isConjugated = True


    def setCarbonylAtoms(self, mol, vRings ):
        for i in range(len(mol.atoms)):
            atom = mol.atoms[i]

            if atom.element.lower() == "c" and atom.num_linkages == 3:
                for j in range(atom.num_linkages):
                    atom_linked = mol.atoms[atom.linkage[j]]

                    if atom.bondOrder[j] == 2 and ((atom_linked.element.lower() == "o") or (atom_linked.element.lower() == "s")):
                        atom.isCarbonyl = True
                        atom_linked.isCarbonyl = True

                        if atom.numCarbonAtoms == 2:
                            atom.isKetone = True
                            atom_linked.isKetone = True

                        if atom_linked.element.lower() == "s":
                            atom_linked.isThiocarbonylS = True

                        break

        hasCarbonyl = False

        for i in range(len(vRings)):
            hasCarbonyl = False
            array = vRings[i].path.split(" ")

            for j in range(len(array)):
                if mol.atoms[int(array[j])].isCarbonyl:
                    hasCarbonyl = True
                    break

            if hasCarbonyl:
                for j in range(len(array)):
                    mol.atoms[int(array[j])].ringHasCabonyl = True


    def setAmide(self, mol, vRings ):
        for i in range(len(mol.atoms)):
            atom = mol.atoms[i]

            if atom.element.lower() == "n":
                for j in range(atom.num_linkages):
                    if atom.bondType[j].lower() == "am":
                        atom.isAmide = True
                        mol.atoms[atom.linkage[j]].isAmide = True

                    if atom.bondOrder[j] == 1 and mol.atoms[atom.linkage[j]].element.lower() == "c":
                        if self.hasDoubleBondedOxygen( mol, atom.linkage[j] ) or self.hasDoubleBondedSulfur( mol, atom.linkage[j]):
                            atom.isAmide = True
                            mol.atoms[atom.linkage[j]].isAmide = True

        hasAmide = False

        for i in range(len(vRings)):
            hasAmide = False
            array = vRings[i].path.split(" ")

            for j in range(len(array)):
                if mol.atoms[int(array[j])].isAmide:
                    hasAmide = True
                    break

            if hasAmide:
                for j in range(len(array)):
                    mol.atoms[int(array[j])].ringHasAmide = True


    def setImineCarbon(self, mol ):
        num_atoms = len(mol.atoms)

        for i in range(num_atoms):
            atom = mol.atoms[i]

            if atom.element.lower() == "c" and atom.num_linkages == 3 and atom.numNitrogenAtoms >= 1:
                for j in range(atom.num_linkages):
                    if atom.bondOrder[j] == 2 or atom.bondType[j].lower() == "ar" and mol.atoms[atom.linkage[j]].element.lower() == "n":
                        atom.isImineCarbon = True
                        break

    def setOxidesLinkedToPhosphorousOrSulfur(self, mol ):
        num_atoms = len(mol.atoms)

        for i in range(num_atoms):
            atom = mol.atoms[i]

            if atom.element.lower() == "p" and atom.numOxygenAtoms >= 1:
                for j in range(atom.num_linkages):
                    atom_linked = mol.atoms[atom.linkage[j]]

                    if atom_linked.element.lower() == "o":
                        atom_linked.isPOxide = True

            elif atom.element.lower() == "s" and atom.numOxygenAtoms >= 1:
                for j in range(atom.num_linkages):
                    atom_linked = mol.atoms[atom.linkage[j]]

                    if atom_linked.element.lower() == "o":
                        atom_linked.isSOxide = True

    def setEster(self, mol ):
        num_atoms = len(mol.atoms)

        for i in range(num_atoms):
            atom = mol.atoms[i]

            if atom.element.lower() == "c" and atom.num_linkages == 3 and atom.numOxygenAtoms >= 2:
                atom.isEster = True

                for j in range(3):
                    atom_linked = mol.atoms[atom.linkage[j]]

                    if atom_linked.element.lower() == "o":
                        atom_linked.isEster = True

    def setAtomProtonationState(self, mol ):
        for i in range(len(mol.atoms)):
            atom = mol.atoms[i]

            if atom.element.lower() == "n" and atom.isAromatic == False and atom.num_linkages == 4:
                atom.isProtonatedNitrogen = True
            elif atom.element.lower() == "n" and atom.isAromatic == False and atom.num_linkages == 3 and atom.numHydrogenAtoms >= 1 and atom.totalBondOrder == 4:
                atom.isProtonatedNitrogen = True
            elif atom.element.lower() == "n" and atom.isAromatic == True and atom.num_linkages == 3 and atom.numHydrogenAtoms == 1:
                atom.isProtonatedNitrogen = True
            elif atom.element.lower() == "s" and atom.num_linkages == 1 and atom.numHydrogenAtoms == 0:
                atom.isDeprotonatedSulfur = True
            elif atom.element.lower() == "o" and atom.num_linkages == 1 and atom.bondOrder[0] == 1 and atom.numHydrogenAtoms == 0:
                atom.isDeprotonatedOxygen = True

            if atom.isProtonatedNitrogen == True:
                for j in range(atom.num_linkages):
                    if mol.atoms[atom.linkage[j]].isImineCarbon:
                        mol.atoms[atom.linkage[j]].isPortonatedImineGroup = True
                        break

            if atom.isDeprotonatedOxygen == True:
                atom_linked = mol.atoms[atom.linkage[0]]

                for j in range(atom_linked.num_linkages):
                    if mol.atoms[atom_linked.linkage[j]].element.lower() == "o" and mol.atoms[atom_linked.linkage[j]].num_linkages == 1 and \
                        mol.atoms[atom_linked.linkage[j]].numHydrogenAtoms == 0:
                        mol.atoms[atom_linked.linkage[j]].isDeprotonatedOxygen = True



    def hasDoubleBondedCarbon(self, mol, atm_index):
            atom = mol.atoms[atm_index]
            i = 0
            while i < atom.num_linkages:
                if atom.bondOrder[i] == 2 and mol.atoms[atom.linkage[i]].element.lower() == "C".lower():
                    return True
                i += 1
            return False

    def hasDoubleBondedOxygen(self, mol, atm_index):
        atom = mol.atoms[atm_index]
        i = 0
        while i < atom.num_linkages:
            if atom.bondOrder[i] == 2 and mol.atoms[atom.linkage[i]].element.lower() == "O".lower():
                return True
            i += 1
        return False

    def hasDoubleBondedSulfur(self, mol, atm_index):
        atom = mol.atoms[atm_index]
        i = 0
        while i < atom.num_linkages:
            if atom.bondOrder[i] == 2 and mol.atoms[atom.linkage[i]].element.lower() == "S".lower():
                return True
            i += 1
        return False

    def getInitialEdgesUsingBonds(self, mol):
        vEdges = []
        edge = Edge()
        i = int()
        j = int()
        b = 0
        while b < len(mol.bonds):
            edge = Edge()
            i = mol.bonds[b].i
            j = mol.bonds[b].j
            if i < j:
                edge.i = i
                edge.j = j
            else:
                edge.i = j
                edge.j = i
            edge.path = str(edge.i) + " " + str(edge.j)
            vEdges.append(edge)
            b += 1
        return vEdges

    ## =====================================================
    ## The cutoff bond lengths are based on
    ## the maximum distances for given atom pairs in CGenFF (par_all36_cgen) and
    ## the Handbook of Chemistry and Physics (63rd edition mostly, some information from the 75th        edition)
    ## =====================================================
    def getCutoffDistance(self, atom_i, atom_j):
        if (atom_i.element.lower() == "C".lower() and atom_j.element.lower() == "CL".lower()) or (atom_i.element.lower() == "CL".lower() and atom_j.element.lower() == "C".lower()):
            return 1.89
        elif (atom_i.element.lower() == "N".lower() and atom_j.element.lower() == "P".lower()) or (atom_i.element.lower() == "P".lower() and atom_j.element.lower() == "N".lower()):
            return 1.89
        elif (atom_i.element.lower() == "C".lower() and atom_j.element.lower() == "S".lower()) or (atom_i.element.lower() == "S".lower() and atom_j.element.lower() == "C".lower()):
            return 1.94
        elif (atom_i.element.lower() == "C".lower() and atom_j.element.lower() == "P".lower()) or (atom_i.element.lower() == "P".lower() and atom_j.element.lower() == "C".lower()):
            return 1.99
        elif (atom_i.element.lower() == "C".lower() and atom_j.element.lower() == "BR".lower()) or (atom_i.element.lower() == "BR".lower() and atom_j.element.lower() == "C".lower()):
            return 2.07
        elif atom_i.element.lower() == "S".lower() and atom_j.element.lower() == "S".lower():
            return 2.13
        elif (atom_i.element.lower() == "C".lower() and atom_j.element.lower() == "I".lower()) or (atom_i.element.lower() == "I".lower() and atom_j.element.lower() == "C".lower()):
            return 2.22
        elif (atom_i.element.lower() == "C".lower() and atom_j.element.lower() == "SI".lower()) or (atom_i.element.lower() == "SI".lower() and atom_j.element.lower() == "C".lower()):
            return 1.96
        elif (atom_i.element.lower() == "C".lower() and atom_j.element.lower() == "SE".lower()) or (atom_i.element.lower() == "SE".lower() and atom_j.element.lower() == "C".lower()):
            return 1.99
        elif (atom_i.element.lower() == "C".lower() and atom_j.element.lower() == "AS".lower()) or (atom_i.element.lower() == "AS".lower() and atom_j.element.lower() == "C".lower()):
            return 2.08
        elif (atom_i.element.lower() == "C".lower() and atom_j.element.lower() == "TE".lower()) or (atom_i.element.lower() == "TE".lower() and atom_j.element.lower() == "C".lower()):
            return 2.15
        else:
            return 1.85

    def getMaximumNumberofCovalentBonds(self, element):
        if element.lower() == "C".lower():
            return 4
        elif element.lower() == "N".lower():
            return 4
        elif element.lower() == "O".lower():
            return 2
        elif element.lower() == "P".lower():
            return 5
        elif element.lower() == "S".lower():
            return 6
        else:
            return 1

    def getNumRemainingVertex(self, vMol):
        num = 0
        num_atoms = len(vMol)
        i = 0
        while i < num_atoms:
            if vMol[i].isAvailable:
                num += 1
            i += 1
        return num

    def chooseVertex(self, vMol):
        num_atoms = len(vMol)
        index = -1
        min_num_linkages = 10
        i = 0
        while i < num_atoms:
            if vMol[i].isAvailable:
                if vMol[i].num_linkages < min_num_linkages:
                    min_num_linkages = vMol[i].num_linkages
                    index = i
            i += 1
        return index

    def chooseSingleLinkedVertex(self, vMol):
        num_atoms = len(vMol)
        index = -1
        i = 0
        while i < num_atoms:
            if vMol[i].isAvailable:
                if vMol[i].num_linkages == 1:
                    return i
            i += 1
        return index

    def removeLinkage(self, vMol, index):
        i = 0
        while i < 5:
            if vMol[index].linkage[i] != -1:
                j = 0
                while j < 5:
                    if vMol[vMol[index].linkage[i]].linkage[j] == index:
                        vMol[vMol[index].linkage[i]].linkage[j] = -1
                        vMol[vMol[index].linkage[i]].num_linkages -= 1
                        break
                    j += 1
            i += 1

    def replaceLinkage(self, vMol, vertex, index_i, index_j):
        i = 0
        while i < 5:
            if vMol[index_i].linkage[i] == vertex:
                vMol[index_i].linkage[i] = index_j
            if vMol[index_j].linkage[i] == vertex:
                vMol[index_j].linkage[i] = index_i
            i += 1

    def getNumLinkages(self, atom):
        num_linkages = 0
        i = 0
        while i < 5:
            if atom.linkage[i] != -1:
                num_linkages += 1
            i += 1
        return num_linkages

    def remove(self, vertex, vMol, vEdges, vRings):
        edge = Edge()
        edge_i = Edge()
        edge_j = Edge()
        isRing = False
        isEdge = False
        vertex_related_edges = ""

        for i in range(len(vEdges)) : 
            if vEdges[i].i == vertex or vEdges[i].j == vertex:
                vertex_related_edges += str(i) + " "

        vertex_related_edges = vertex_related_edges.strip()
        array = vertex_related_edges.split(" ")

        if len(vertex_related_edges) != 0 and len(array) == 1:
            self.removeLinkage(vMol, vertex)

        if len(vertex_related_edges) != 0 and len(array) >=2:
            for i in range(len(array)):
                edge_i = vEdges[int(array[i])]

                for j in range(i+1, len(array)):
                    edge_j = vEdges[int(array[j])]

                    isRing = False
                    isEdge = False

                    if edge_i.i == edge_j.i and edge_i.j == edge_j.j:
                        isRing = True

                        edge = Edge()
                        edge.i = vertex
                        edge.j = vertex

                    elif edge_i.i == edge_j.i:
                        isEdge = True
                        edge = Edge()

                        if edge_i.j < edge_j.j :
                            edge.i = edge_i.j
                            edge.j = edge_j.j
                        else:
                            edge.i = edge_j.j
                            edge.j = edge_i.j
                        self.replaceLinkage(vMol, vertex, edge_i.j, edge_j.j)

                    elif edge_i.i == edge_j.j:
                        isEdge = True
                        edge = Edge()
                        
                        if edge_i.j < edge_j.i :
                            edge.i = edge_i.j
                            edge.j = edge_j.i
                        else:
                            edge.i = edge_j.i
                            edge.j = edge_i.j
                        
                        self.replaceLinkage( vMol, vertex, edge_i.j, edge_j.i )

                    elif edge_i.j == edge_j.i:
                        isEdge = True
                        edge = Edge()
                        
                        if edge_i.i < edge_j.j:
                            edge.i = edge_i.i
                            edge.j = edge_j.j
                        
                        else:
                            edge.i = edge_j.j
                            edge.j = edge_i.i
                        
                        self.replaceLinkage( vMol, vertex, edge_i.i, edge_j.j )
                    
                    elif edge_i.j == edge_j.j:
                        isEdge = True
                        edge = Edge()
                        
                        if edge_i.i < edge_j.i :
                            edge.i = edge_i.i
                            edge.j = edge_j.i
    
                        else:
                            edge.i = edge_j.i
                            edge.j = edge_i.i
                        
                        self.replaceLinkage( vMol, vertex, edge_i.i, edge_j.i )

                    if isRing == True:
                        edge.path = self.removeCommonVertex( edge_i.path, edge_j.path ) + " " + edge_j.path
                        edge.path = edge.path.strip()
                        vRings.append( edge )

                    
                    if isEdge == True:
                        edge.path = self.removeCommonVertex( edge_i.path, edge_j.path ) + " " + edge_j.path
                        edge.path = edge.path.strip()
                        vEdges.append( edge )

        if  len(vertex_related_edges) != 0 and  len(array) >= 1:
            for i in range(len(array)):
                vEdges[int( array[i] ) ].path = None
            
            for i in range(len(vEdges)):
                if vEdges[i].path == None:
                    del vEdges[i]
                    i -= 1
                
        self.removeLinkage(vMol, vertex)

        vMol[vertex].isAvailable = False
                        

    def removeCommonVertex(self, path_i, path_j):
        array_i = path_i.split(" ")
        array_j = path_j.split(" ")
        new_path = ""
        vertex_i = int()
        vertex_j = int()
        isCommonVertex = False
        i = 0
        while len(array_i):
            vertex_i = int(array_i[i])
            isCommonVertex = False
            while len(array_j):
                vertex_j = int(array_j[j])
                if vertex_i == vertex_j:
                    isCommonVertex = True
                    break
                j += 1
            if not isCommonVertex:
                new_path += array_i[i] + " "
            i += 1
        return new_path.strip()

    def mergePaths(self, path_i, path_j):
        array_i = path_i.split(" ")
        array_j = path_j.split(" ")
        new_path = ""
        vertex_i = int()
        vertex_j = int()
        isCommonVertex = False
        i = 0
        while len(array_i):
            vertex_i = int(array_i[i])
            isCommonVertex = False
            while len(array_j):
                vertex_j = int(array_j[j])
                if vertex_i == vertex_j:
                    isCommonVertex = True
                    break
                j += 1
            if not isCommonVertex:
                new_path += array_i[i] + " "
            i += 1
        if 0 != len(length):
            return new_path.strip() + " " + path_j
        else:
            return path_j

    def removeLargeRings(self, vRings, num_atoms):
        array = []
        i = 0
        while i < len(vRings):
            array = vRings[i].path.split(" ")
            if len(array):
                del vRings[i]
                i -= 1
            i += 1

    def clearEdges(self, vEdges, num_edges):
        array = []
        i = 0
        while i < len(vEdges):
            array = vEdges[i].path.split(" ")
            if len(array):
                del vEdges[i]
                i -= 1
            i += 1

    def isRing(self, vMol, ring):
        array = ring.path.split(" ")
        num_atoms = int()
        index_i = int()
        index_j = int()
        dist = float()
        num_linkages = 0
        i = 0
        while i < num_atoms:
            num_linkages = 0
            index_i = int(array[i])
            while j < num_atoms:
                index_j = int(array[j])
                if i != j:
                    dist = getDistance(vMol[index_i].R, vMol[index_j].R)
                    if dist < 1.85:
                        num_linkages += 1
                j += 1
            if num_linkages == 1:
                return False
            i += 1
        return True

    def removeNonRings(self, vMol, vRings):
        i = 0
        while i < len(vRings):
            if not isRing(vMol, vRings[i]):
                del vRings[i]
                i -= 1
            i += 1

    def removeRedundantRings(self, vRings ):

        num_identical_atoms = 0

        for i in range(len(vRings)):
            edge_i = vRings[i]
            array_1 = edge_i.path.split(" ")

            for j in range(i+1, len(vRings)):
                edge_j = vRings[j]
                array_2 = edge_j.path.split(" ")

                if len(array_1) == len(array_2):
                    num_identical_atoms = 0

                    for s1 in range(len(array_1)):
                        for s2 in range(len(array_2)):
                            if int(array_1[s1]) == int(array_2[s2]):
                                num_identical_atoms = num_identical_atoms + 1
                                break

                    if num_identical_atoms == len(array_2):
                        vRings.pop(j)
                        j = j - 1


    def keepSmallestRings(self, vRings ):
        num_identical_atoms = 0
        option = 0
        SubRingIsRemoved = True

        while(SubRingIsRemoved):
            SubRingIsRemoved = False

            for i in range(len(vRings)):
                edge_i = vRings[i]
                array_1 = edge_i.path.split(" ")

                for j in range(i+1, len(vRings)):
                    edge_j = vRings[j]
                    array_2 = edge_j.path.split(" ")

                    if len(array_1) <= len(array_2):
                        path_1 = edge_i.path
                        path_2 = edge_j.path
                        option = 1
                    else:
                        path_1 = edge_j.path
                        path_2 = edge_i.path
                        option = 2

                    array_1 = path_1.split(" ")
                    array_2 = path_2.split(" ")

                    num_identical_atoms = 0

                    for s1 in range(len(array_1)):
                        for s2 in range(len(array_2)):
                            if int(array_1[s1]) == int(array_2[s2]):
                                num_identical_atoms = num_identical_atoms + 1
                                break

                    if num_identical_atoms == len(array_1):
                        SubRingIsRemoved = True
                        if option == 1:
                            vRings.pop(j)
                        elif option == 2:
                            vRings.pop(i)
                        break

                if SubRingIsRemoved:
                    break



    def mergeRings(self, vRings):
        edge_i = Edge()
        edge_j = Edge()
        path_1 = str()
        path_2 = str()
        array_1 = []
        array_2 = []
        mergeRings = True
        while mergeRings:
            mergeRings = False
            while i < len(vRings):
                edge_i = vRings[i]
                path_1 = edge_i.path
                array_1 = path_1.split(" ")
                while j < len(vRings):
                    edge_j = vRings[j]
                    path_2 = edge_j.path
                    array_2 = path_2.split(" ")
                    while len(array_1):
                        while len(array_2):
                            if int(array_1[s1]) == int(array_2[s2]):
                                mergeRings = True
                                vRings[i].path = mergePaths(path_1, path_2)
                                del vRings[j]
                                break
                            s2 += 1
                        if mergeRings:
                            break
                        s1 += 1
                    if mergeRings:
                        break
                    j += 1
                if mergeRings:
                    break
                i += 1





