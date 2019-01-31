#!/usr/bin/env python

import math
from SmallMolecule import SmallMolecule

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

        with open (file_name) as bufferedreader:
            for line in bufferedreader.readlines():
                #  *** read atom information ***
                if line.startsWith("@<TRIPOS>ATOM"):
                    if line.startsWith("@<TRIPOS>BOND"):
                        break
                    #  atom_id / atom_name / x / y / z / atom_type / subst_id / subst_name / charge / status_bit
                    atom = Atom()
                    atom.line = line
                    atom.atom_name = line.strip().split()[1]
                    atom.R[0] = float(line.strip().split()[2])
                    atom.R[1] = float(line.strip().split()[3])
                    atom.R[2] = float(line.strip().split()[4])
                    element = line.strio().split()[5]
                    index = element.indexOf('.')
                    if index != -1:
                        atom.element = element.substring(0, index)
                    else:
                        atom.element = element
                    atom.element = atom.element.upper()
                    mol.atoms.append(atom)

                #  *** read bond information ***
                if line.startsWith("@<TRIPOS>BOND"):
                    if line.startsWith("@<TRIPOS>SUBSTRUCTURE") or (0 == len(line)):
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
                    mol.bonds.add(bond)
        return mol

    def readCGenFFBondParameters(self, file_name):

    def readCGenFFAngleParameters(self, file_name):

    def readCGenFFDihedralParameters(self, file_name):

    def readCGenFFImproperParameters(self, file_name):

    def readAtomPriority(self, file_name):

    def readAtomTypes(self, file_name):

    def setElement(self, mol, hAtomType):

    def readChargeIncrementBond(self, file_name):

    def orderAtoms(self, hAtomPriority, atoms_combi, kind):

    def orderAtomsBonds(self, hAtomPriority, mol):

    def orderAtomsAngles(self, hAtomPriority, mol):

    def setFormalCharge(self, mol):
        num_atoms = len(mol.atoms)
        atom = Atom()
        i = 0
        while i < num_atoms:
            atom = mol.atoms.elementAt(i)
            #  -1
            #  heavy atoms linked to OG2D2, OG312, SG302, OG2P1
            #  PG1 linked to OG2P1
            if atom.atomType.lower() == "OG2D2".lower() or atom.atomType.lower() == "OG312".lower() or atom.atomType.lower() == "SG302".lower() or atom.atomType.lower() == "OG2P1".lower():
                while j < 5:
                    if atom.linkage[j] != -1:
                        if mol.atoms.elementAt(atom.linkage[j]).element.charAt(0) == 'C' or mol.atoms.elementAt(atom.linkage[j]).element.charAt(0) == 'S':
                            mol.atoms.elementAt(atom.linkage[j]).formal_charge = -1
                    j += 1
            if atom.atomType.lower() == "OG2P1".lower():
                if mol.atoms.elementAt(atom.linkage[0]).atomType.lower() == "PG1".lower():
                    #  -1
                    #  PG1 linked to OG2P1
                    mol.atoms.elementAt(atom.linkage[0]).formal_charge = -1
                elif mol.atoms.elementAt(atom.linkage[0]).atomType.lower() == "PG2".lower():
                    #  -2
                    #  PG2 linked to to OG2P1
                    mol.atoms.elementAt(atom.linkage[0]).formal_charge = -2
            #  1
            #  heavy atoms linked to NG2P1, NG3P3, NG2R52, NG3P2, NG3P0, NG2R61, NG3P1
            if atom.atomType.lower() == "NG2P1".lower() or atom.atomType.lower() == "NG3P3".lower() or atom.atomType.lower() == "NG2R52".lower() or atom.atomType.lower() == "NG3P2".lower() or atom.atomType.lower() == "NG3P0".lower() or atom.atomType.lower() == "NG2R61".lower() or atom.atomType.lower() == "NG3P1".lower():
                while j < 5:
                    if atom.linkage[j] != -1:
                        if mol.atoms.elementAt(atom.linkage[j]).element.charAt(0) == 'C':
                            mol.atoms.elementAt(atom.linkage[j]).formal_charge = 1
                    j += 1
            i += 1

    def removeZeroIncrementBond(self, mol):
        #  The charge increment for a bond between two identical atom types is zero by definition.
        i = 0
        while i < len(mol.bonds):
            if mol.atoms.elementAt(mol.bonds.elementAt(i).i).atomType.lower() == mol.atoms.elementAt(mol.bonds.elementAt(i.lower().j).atomType):
                mol.bonds.removeElementAt(i)
                i -= 1
            elif mol.atoms.elementAt(mol.bonds.elementAt(i).i).atomType.lower() == "CG2D1O".lower() and mol.atoms.elementAt(mol.bonds.elementAt(i).j).atomType.lower() == "CG2D2O".lower():
                mol.bonds.removeElementAt(i)
                i -= 1
            elif mol.atoms.elementAt(mol.bonds.elementAt(i).i).atomType.lower() == "CG2D2O".lower() and mol.atoms.elementAt(mol.bonds.elementAt(i).j).atomType.lower() == "CG2D1O".lower():
                mol.bonds.removeElementAt(i)
                i -= 1
            elif mol.atoms.elementAt(mol.bonds.elementAt(i).i).atomType.lower() == "CG2DC1".lower() and mol.atoms.elementAt(mol.bonds.elementAt(i).j).atomType.lower() == "CG2DC2".lower():
                mol.bonds.removeElementAt(i)
                i -= 1
            elif mol.atoms.elementAt(mol.bonds.elementAt(i).i).atomType.lower() == "CG2DC2".lower() and mol.atoms.elementAt(mol.bonds.elementAt(i).j).atomType.lower() == "CG2DC1".lower():
                mol.bonds.removeElementAt(i)
                i -= 1
            elif mol.atoms.elementAt(mol.bonds.elementAt(i).i).atomType.lower() == "CG25C1".lower() and mol.atoms.elementAt(mol.bonds.elementAt(i).j).atomType.lower() == "CG25C2".lower():
                mol.bonds.removeElementAt(i)
                i -= 1
            elif mol.atoms.elementAt(mol.bonds.elementAt(i).i).atomType.lower() == "CG25C2".lower() and mol.atoms.elementAt(mol.bonds.elementAt(i).j).atomType.lower() == "CG25C1".lower():
                mol.bonds.removeElementAt(i)
                i -= 1
            elif mol.atoms.elementAt(mol.bonds.elementAt(i).i).atomType.lower() == "CG251O".lower() and mol.atoms.elementAt(mol.bonds.elementAt(i).j).atomType.lower() == "CG252O".lower():
                mol.bonds.removeElementAt(i)
                i -= 1
            elif mol.atoms.elementAt(mol.bonds.elementAt(i).i).atomType.lower() == "CG252O".lower() and mol.atoms.elementAt(mol.bonds.elementAt(i).j).atomType.lower() == "CG251O".lower():
                mol.bonds.removeElementAt(i)
                i -= 1
            i += 1

    def setChargeIncrementBond(self, mol, hChargeIncrement):
        value = None
        atom_combi = ""
        i = 0
        while i < len(mol.bonds):
            atom_combi = mol.atoms.elementAt(mol.bonds.elementAt(i).i).atomType + " " + mol.atoms.elementAt(mol.bonds.elementAt(i).j).atomType
            value = hChargeIncrement.get(atom_combi)
            if value != None:
                mol.bonds.elementAt(i).increment = Double.parseDouble(value)
            else:
                print "No available bond charge increment: " + atom_combi
            i += 1

    def assignAtomChargesUsingBondIncrement(self, mol):
        a = 0
        while a < len(mol.atoms):
            mol.atoms.elementAt(a).charge_calc = mol.atoms.elementAt(a).formal_charge
            a += 1
        b = 0
        while b < len(mol.bonds):
            mol.atoms.elementAt(mol.bonds.elementAt(b).i).charge_calc -= mol.bonds.elementAt(b).increment
            mol.atoms.elementAt(mol.bonds.elementAt(b).j).charge_calc += mol.bonds.elementAt(b).increment
            b += 1

    #  =============================================================================
    #  Algorithm is based on Hanser et al. J. Chem. Inf. Comput. Sci. 1996, 36, 1146
    #  A new algorithm for exhaustive ring perception in a molecular graph
    #  =============================================================================
    def detectRings(self, mol):
        atom = Atom()
        edge = Edge()
        frag = Fragment()
        vEdges = Vector()
        vRings = Vector()
        #  *** generate P-Graph ***
        vEdges = getInitialEdgesUsingBonds(mol)
        while getNumRemainingVertex(mol.atoms) != 0:
            clearEdges(vEdges, 7)
            remove(vertex, mol.atoms, vEdges, vRings)
        removeNonRings(mol.atoms, vRings)
        removeRedundantRings(vRings)
        keepSmallestRings(vRings)
        return vRings

    def setRingAtoms(self, mol, vRings):
        array = []
        index = int()
        i = 0
        while i < len(vRings):
            array = vRings.elementAt(i).path.split(" ")
            while len(array):
                index = Integer.parseInt(array[j])
                mol.atoms.elementAt(index).isRingAtom = True
                while k < 3:
                    if mol.atoms.elementAt(index).ringIndex[k] == -1:
                        mol.atoms.elementAt(index).ringIndex[k] = i
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
            array = vRings.elementAt(i).path.split(" ")
            # print  vRings.elementAt(i).path ;
            #  *** check planarity ***
            if len(array):
                if isPlanarRing(mol, array):
                    aRingAromacity[i] = True
                    while len(array):
                        index = Integer.parseInt(array[j])
                        mol.atoms.elementAt(index).isAromatic = True
                        if mol.atoms.elementAt(index).isCarbonyl:
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
                index = Integer.parseInt(array[j])
                mol.atoms.elementAt(index).isPyrrole = is_pyrrole
                mol.atoms.elementAt(index).isNeutralPyridine = is_neutral_pyridine
                mol.atoms.elementAt(index).isProtonatedPyridine = is_protonated_pyridine
                mol.atoms.elementAt(index).isProtonatedPyrimidine = is_protonated_pyrimidine
                mol.atoms.elementAt(index).isFuran = is_furan
                mol.atoms.elementAt(index).isFuranose = is_furanose
                j += 1
            i += 1

    def isPlanarRing(self, mol, array):
        atom = Atom()
        i = 0
        while len(array):
            atom = mol.atoms.elementAt(Integer.parseInt(array[i]))
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
                atom = mol.atoms.elementAt(Integer.parseInt(array[i]))
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
                atom = mol.atoms.elementAt(Integer.parseInt(array[i]))
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
                atom = mol.atoms.elementAt(Integer.parseInt(array[i]))
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
                atom = mol.atoms.elementAt(Integer.parseInt(array[i]))
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
                atom = mol.atoms.elementAt(Integer.parseInt(array[i]))
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
                index = Integer.parseInt(array[i])
                atom = mol.atoms.elementAt(index)
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
            while j < 5:
                if mol.atoms.elementAt(i).linkage[j] != -1:
                    if mol.atoms.elementAt(mol.atoms.elementAt(i).linkage[j]).element.lower() == "H".lower():
                        num_hydrogen_atoms += 1
                j += 1
            mol.atoms.elementAt(i).numHydrogenAtoms = num_hydrogen_atoms
            i += 1

    def setNumCarbonAtoms(self, mol):
        num_carbon_atoms = 0
        i = 0
        while i < len(mol.atoms):
            num_carbon_atoms = 0
            while j < 5:
                if mol.atoms.elementAt(i).linkage[j] != -1:
                    if mol.atoms.elementAt(mol.atoms.elementAt(i).linkage[j]).element.lower() == "C".lower():
                        num_carbon_atoms += 1
                j += 1
            mol.atoms.elementAt(i).numCarbonAtoms = num_carbon_atoms
            i += 1

    def setNumNitrogenAtoms(self, mol):
        num_nigrogen_atoms = 0
        i = 0
        while i < len(mol.atoms):
            num_nigrogen_atoms = 0
            while j < 5:
                if mol.atoms.elementAt(i).linkage[j] != -1:
                    if mol.atoms.elementAt(mol.atoms.elementAt(i).linkage[j]).element.lower() == "N".lower():
                        num_nigrogen_atoms += 1
                j += 1
            mol.atoms.elementAt(i).numNitrogenAtoms = num_nigrogen_atoms
            i += 1

    def setNumOxygenAtoms(self, mol):
        num_oxygen_atoms = 0
        i = 0
        while i < len(mol.atoms):
            num_oxygen_atoms = 0
            while j < 5:
                if mol.atoms.elementAt(i).linkage[j] != -1:
                    if mol.atoms.elementAt(mol.atoms.elementAt(i).linkage[j]).element.lower() == "O".lower():
                        num_oxygen_atoms += 1
                j += 1
            mol.atoms.elementAt(i).numOxygenAtoms = num_oxygen_atoms
            i += 1

    def setNumSulfurAtoms(self, mol):
        num_sulfur_atoms = 0
        i = 0
        while i < len(mol.atoms):
            num_sulfur_atoms = 0
            while j < 5:
                if mol.atoms.elementAt(i).linkage[j] != -1:
                    if mol.atoms.elementAt(mol.atoms.elementAt(i).linkage[j]).element.lower() == "S".lower():
                        num_sulfur_atoms += 1
                j += 1
            mol.atoms.elementAt(i).numSulfurAtoms = num_sulfur_atoms
            i += 1

    def setMethylAtoms(self, mol):
        i = 0
        while i < len(mol.atoms):
            if mol.atoms.elementAt(i).element.lower() == "C".lower() and mol.atoms.elementAt(i).numHydrogenAtoms == 3:
                mol.atoms.elementAt(i).isMethyl = True
            i += 1

    def isAtomInFiveMemberedRing(self, mol, atm_index):
        atom = mol.atoms.elementAt(atm_index)
        i = 0
        while i < 3:
            if atom.numRingAtoms[i] == 5:
                return True
            i += 1
        return False

    def isAtomInSixMemberedRing(self, mol, atm_index):
        atom = mol.atoms.elementAt(atm_index)
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
            while j < 5:
                if mol.atoms.elementAt(mol.bonds.elementAt(i).i).linkage[j] == mol.bonds.elementAt(i).j:
                    if mol.bonds.elementAt(i).bondType.charAt(0) == '1':
                        mol.atoms.elementAt(mol.bonds.elementAt(i).i).bondOrder[j] = 1
                    elif mol.bonds.elementAt(i).bondType.charAt(0) == '2':
                        mol.atoms.elementAt(mol.bonds.elementAt(i).i).bondOrder[j] = 2
                    elif mol.bonds.elementAt(i).bondType.charAt(0) == '3':
                        mol.atoms.elementAt(mol.bonds.elementAt(i).i).bondOrder[j] = 3
                    elif mol.bonds.elementAt(i).bondType.lower() == "am".lower():
                        mol.atoms.elementAt(mol.bonds.elementAt(i).i).bondOrder[j] = 1
                    break
                j += 1
            while j < 5:
                if mol.atoms.elementAt(mol.bonds.elementAt(i).j).linkage[j] == mol.bonds.elementAt(i).i:
                    if mol.bonds.elementAt(i).bondType.charAt(0) == '1':
                        mol.atoms.elementAt(mol.bonds.elementAt(i).j).bondOrder[j] = 1
                    elif mol.bonds.elementAt(i).bondType.charAt(0) == '2':
                        mol.atoms.elementAt(mol.bonds.elementAt(i).j).bondOrder[j] = 2
                    elif mol.bonds.elementAt(i).bondType.charAt(0) == '3':
                        mol.atoms.elementAt(mol.bonds.elementAt(i).j).bondOrder[j] = 3
                    elif mol.bonds.elementAt(i).bondType.lower() == "am".lower():
                        mol.atoms.elementAt(mol.bonds.elementAt(i).j).bondOrder[j] = 1
                    break
                j += 1
            i += 1
        i = 0
        while i < len(mol.atoms):
            atom = mol.atoms.elementAt(i)
            if atom.isAromatic and atom.element.lower() == "C".lower() and hasDoubleBondedOxygen(mol, i) and atom.num_linkages == 3:
                while j < 3:
                    if not mol.atoms.elementAt(atom.linkage[j]).element.lower() == "O".lower():
                        atom.bondOrder[j] = 1
                        while k < 5:
                            atom_linked = mol.atoms.elementAt(atom.linkage[j])
                            if atom_linked.linkage[k] == i:
                                atom_linked.bondOrder[k] = 1
                                break
                            k += 1
                    j += 1
            i += 1
        i = 0
        while i < len(mol.atoms):
            atom = mol.atoms.elementAt(i)
            if atom.isAromatic and atom.element.lower() == "C".lower() and atom.num_linkages == 3 and getNumUnassignedBondOrders(atom) == 1:
                while j < 3:
                    if atom.bondOrder[j] == 0:
                        atom.bondOrder[j] = 2
                        while k < 5:
                            atom_linked = mol.atoms.elementAt(atom.linkage[j])
                            if atom_linked.linkage[k] == i:
                                atom_linked.bondOrder[k] = 2
                                break
                            k += 1
                        break
                    j += 1
            i += 1
        i = 0
        while i < len(mol.atoms):
            atom = mol.atoms.elementAt(i)
            if atom.isAromatic and atom.element.lower() == "C".lower() and atom.num_linkages == 3 and getNumUnassignedBondOrders(atom) == 2:
                #  double bond for one bond
                while j < 3:
                    if atom.bondOrder[j] == 0:
                        atom.bondOrder[j] = 2
                        while k < 5:
                            atom_linked = mol.atoms.elementAt(atom.linkage[j])
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
                            atom_linked = mol.atoms.elementAt(atom.linkage[j])
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
        while i < num_atoms:
            atom = mol.atoms.elementAt(i)
            total_bond_order = 0
            while j < 5:
                total_bond_order += atom.bondOrder[j]
                j += 1
            mol.atoms.elementAt(i).totalBondOrder = total_bond_order
            i += 1

    def setAngleList(self, mol):
        atmgrp = LinkedAtomGroup()
        array = []
        i = 0
        while i < len(mol.bonds):
            while j < len(mol.bonds):
                if mol.bonds.elementAt(i).i == mol.bonds.elementAt(j).j:
                    atmgrp = LinkedAtomGroup()
                    atmgrp.i = mol.bonds.elementAt(j).i
                    atmgrp.j = mol.bonds.elementAt(i).j
                    atmgrp.linked_atoms = "" + mol.bonds.elementAt(j).i + "-" + mol.bonds.elementAt(i).i + "-" + mol.bonds.elementAt(i).j
                    atmgrp.linked_bonds = "" + j + "-" + i
                    mol.angles.add(atmgrp)
                elif mol.bonds.elementAt(i).j == mol.bonds.elementAt(j).i:
                    atmgrp = LinkedAtomGroup()
                    atmgrp.i = mol.bonds.elementAt(i).i
                    atmgrp.j = mol.bonds.elementAt(j).j
                    atmgrp.linked_atoms = "" + mol.bonds.elementAt(i).i + "-" + mol.bonds.elementAt(i).j + "-" + mol.bonds.elementAt(j).j
                    atmgrp.linked_bonds = "" + i + "-" + j
                    mol.angles.add(atmgrp)
                elif mol.bonds.elementAt(i).i == mol.bonds.elementAt(j).i:
                    atmgrp = LinkedAtomGroup()
                    atmgrp.i = mol.bonds.elementAt(j).j
                    atmgrp.j = mol.bonds.elementAt(i).j
                    atmgrp.linked_atoms = "" + mol.bonds.elementAt(j).j + "-" + mol.bonds.elementAt(i).i + "-" + mol.bonds.elementAt(i).j
                    atmgrp.linked_bonds = "" + j + "-" + i
                    mol.angles.add(atmgrp)
                elif mol.bonds.elementAt(i).j == mol.bonds.elementAt(j).j:
                    atmgrp = LinkedAtomGroup()
                    atmgrp.i = mol.bonds.elementAt(i).i
                    atmgrp.j = mol.bonds.elementAt(j).i
                    atmgrp.linked_atoms = "" + mol.bonds.elementAt(i).i + "-" + mol.bonds.elementAt(i).j + "-" + mol.bonds.elementAt(j).i
                    atmgrp.linked_bonds = "" + i + "-" + j
                    mol.angles.add(atmgrp)
                j += 1
            i += 1
        #
        # 		print  len(mol.angles) ;
        # 		for( int i = 0 ; i < len(mol.angles) ; i++ ){
        # 			for( int j = i+1 ; j < len(mol.angles) ; j++ ){
        # 				array = mol.angles.elementAt(j).linked_atoms.split( "-" );
        # 				if( mol.angles.elementAt(i).linked_atoms.lower() ==  mol.angles.elementAt(j.lower().linked_atoms ) ){
        # 					mol.angles.removeElementAt(j);
        # 					j--;
        # 				}
        # 			}
        # 		}
        # 		print  len(vAngle) ;
        #

    def setDihderalList(self, mol):
        atmgrp = LinkedAtomGroup()
        array = []
        i = 0
        while i < len(mol.angles):
            array = mol.angles.elementAt(i).linked_bonds.split("-")
            while j < len(mol.bonds):
                if j != Integer.parseInt(array[0]) and j != Integer.parseInt(array[1]):
                    if mol.angles.elementAt(i).i == mol.bonds.elementAt(j).j:
                        atmgrp = LinkedAtomGroup()
                        atmgrp.i = mol.bonds.elementAt(j).i
                        atmgrp.j = mol.angles.elementAt(i).j
                        atmgrp.linked_atoms = mol.bonds.elementAt(j).i + "-" + mol.angles.elementAt(i).linked_atoms
                        atmgrp.linked_bonds = j + "-" + mol.angles.elementAt(i).linked_bonds
                        mol.dihedrals.add(atmgrp)
                    elif mol.angles.elementAt(i).j == mol.bonds.elementAt(j).i:
                        atmgrp = LinkedAtomGroup()
                        atmgrp.i = mol.angles.elementAt(i).i
                        atmgrp.j = mol.bonds.elementAt(j).j
                        atmgrp.linked_atoms = mol.angles.elementAt(i).linked_atoms + "-" + mol.bonds.elementAt(j).j
                        atmgrp.linked_bonds = mol.angles.elementAt(i).linked_bonds + "-" + j
                        mol.dihedrals.add(atmgrp)
                    elif mol.angles.elementAt(i).i == mol.bonds.elementAt(j).i:
                        atmgrp = LinkedAtomGroup()
                        atmgrp.i = mol.bonds.elementAt(j).j
                        atmgrp.j = mol.angles.elementAt(i).j
                        atmgrp.linked_atoms = mol.bonds.elementAt(j).j + "-" + mol.angles.elementAt(i).linked_atoms
                        atmgrp.linked_bonds = j + "-" + mol.angles.elementAt(i).linked_bonds
                        mol.dihedrals.add(atmgrp)
                    elif mol.angles.elementAt(i).j == mol.bonds.elementAt(j).j:
                        atmgrp = LinkedAtomGroup()
                        atmgrp.i = mol.angles.elementAt(i).i
                        atmgrp.j = mol.bonds.elementAt(j).i
                        atmgrp.linked_atoms = mol.angles.elementAt(i).linked_atoms + "-" + mol.bonds.elementAt(j).i
                        atmgrp.linked_bonds = mol.angles.elementAt(i).linked_bonds + "-" + j
                        mol.dihedrals.add(atmgrp)
                j += 1
            i += 1
        # print  len(vDihedral) ;
        i = 0
        while i < len(mol.dihedrals):
            while j < len(mol.dihedrals):
                if mol.dihedrals.elementAt(i).linked_atoms.lower() == mol.dihedrals.elementAt(j.lower().linked_atoms):
                    mol.dihedrals.removeElementAt(j)
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
            atom = mol.atoms.elementAt(i)
            atom_type = atom.atomType
            if atom.element.lower() == "C".lower() and atom.num_linkages == 3:
                if not atom.isRingAtom and getNumberOfLinkedRingAtoms(mol, i) <= 1:
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
                elif atom.isRingAtom and getNumberOfLinkedRingAtoms(mol, i) == 2:
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
                if not atom.isRingAtom and getNumberOfLinkedRingAtoms(mol, i) <= 1:
                    if hasDoubleBondedOxygen(mol, i) and atom.numCarbonAtoms == 1 and atom.numOxygenAtoms == 2:
                        add_to_list = True
                    elif getNumberOfLinkedRingAtoms(mol, i) == 1 and atom.numCarbonAtoms == 1 and atom.numHydrogenAtoms == 2:
                        add_to_list = True
            if add_to_list:
                atmgrp = LinkedAtomGroup()
                atmgrp.linked_atoms = "" + i
                while j < 3:
                    atmgrp.linked_atoms += "-" + atom.linkage[j]
                    j += 1
                mol.impropers.add(atmgrp)
            i += 1


    def getNumberOfLinkedRingAtoms(self, mol, atm_index):
        atom = mol.atoms.elementAt(atm_index)
        num_ring_atoms = 0
        i = 0
        while i < atom.num_linkages:
            if mol.atoms.elementAt(atom.linkage[i]).isRingAtom:
                num_ring_atoms += 1
            i += 1
        return num_ring_atoms

    def getFirstExocylicAtom(self, mol, atm_index):
        atom = mol.atoms.elementAt(atm_index)
        i = 0
        while i < atom.num_linkages:
            if not mol.atoms.elementAt(atom.linkage[i]).isRingAtom:
                return mol.atoms.elementAt(atom.linkage[i]).element
            i += 1
        return "X"

    def sortAtoms(self, atm_grp, kind):
        #  ORDER OF PREFERENCE FOR SORTING PARAMETERS:
        #  C < N < O < P < S < HALOGENS (LOW TO HIGH Z) < MISC. (BY Z) < H
        #  ATOMS TYPES WITHIN THE SAME ELEMENT ARE SORTED ALPHABETICALLY
        atoms = ""
        if kind == 1:
            #  BOND
            #  THE LOWEST PRIORITY ATOM ALWAYS COMES FIRST
        elif kind == 2:
            #  ANGLE
        elif kind == 3:
            #  DIHEDRAL
        elif kind == 4:
            #  IMPROPER
        return atoms

    def setConjugatedAtoms(self, mol):

    public void setCarbonylAtoms( SmallMolecule mol, Vector<Edge> vRings ){
    public void setAmide( SmallMolecule mol, Vector<Edge> vRings ){
    public void setImineCarbon( SmallMolecule mol ){
    public void setOxidesLinkedToPhosphorousOrSulfur( SmallMolecule mol ){
    public void setEster( SmallMolecule mol ){
    public void setAtomProtonationState( SmallMolecule mol ){
    public boolean hasDoubleBondedCarbon( SmallMolecule mol, int atm_index ){
    public boolean hasDoubleBondedOxygen( SmallMolecule mol, int atm_index ){
    public boolean hasDoubleBondedSulfur( SmallMolecule mol, int atm_index ){
    public Vector<Edge> getInitialEdgesUsingBonds( SmallMolecule mol ){
    public double getCutoffDistance( Atom atom_i, Atom atom_j ){
    public int getMaximumNumberofCovalentBonds( String element ){
    public int getNumRemainingVertex( Vector<Atom> vMol ){
    public int chooseVertex( Vector<Atom> vMol ){
    public int chooseSingleLinkedVertex( Vector<Atom> vMol ){
    public void removeLinkage( Vector<Atom> vMol, int index ){
    public void replaceLinkage( Vector<Atom> vMol, int vertex, int index_i, int index_j ){
    public int getNumLinkages( Atom atom ){
    public void remove( int vertex, Vector<Atom> vMol, Vector<Edge> vEdges, Vector<Edge> vRings ){
    public String removeCommonVertex( String path_i, String path_j ){
    public String mergePaths( String path_i, String path_j ){
    public void removeLargeRings( Vector<Edge> vRings, int num_atoms ){
    public void clearEdges( Vector<Edge> vEdges, int num_edges ){
    public boolean isRing( Vector<Atom> vMol, Edge ring ){
    public void removeNonRings( Vector<Atom> vMol, Vector<Edge> vRings ){
    public void removeRedundantRings( Vector<Edge> vRings ){
    public void keepSmallestRings( Vector<Edge> vRings ){
    public void mergeRings( Vector<Edge> vRings ){





