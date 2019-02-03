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



    def hasDoubleBondedCarbon(self, mol, atm_index):
        atom = mol.atoms.elementAt(atm_index)
        i = 0
        while i < atom.num_linkages:
            if atom.bondOrder[i] == 2 and mol.atoms.elementAt(atom.linkage[i]).element.lower() == "C".lower():
                return True
            i += 1
        return False

    def hasDoubleBondedOxygen(self, mol, atm_index):
        atom = mol.atoms.elementAt(atm_index)
        i = 0
        while i < atom.num_linkages:
            if atom.bondOrder[i] == 2 and mol.atoms.elementAt(atom.linkage[i]).element.lower() == "O".lower():
                return True
            i += 1
        return False

    def hasDoubleBondedSulfur(self, mol, atm_index):
        atom = mol.atoms.elementAt(atm_index)
        i = 0
        while i < atom.num_linkages:
            if atom.bondOrder[i] == 2 and mol.atoms.elementAt(atom.linkage[i]).element.lower() == "S".lower():
                return True
            i += 1
        return False

    def getInitialEdgesUsingBonds(self, mol):
        vEdges = Vector()
        edge = Edge()
        i = int()
        j = int()
        b = 0
        while b < len(mol.bonds):
            edge = Edge()
            i = mol.bonds.elementAt(b).i
            j = mol.bonds.elementAt(b).j
            if i < j:
                edge.i = i
                edge.j = j
            else:
                edge.i = j
                edge.j = i
            edge.path = edge.i + " " + edge.j
            vEdges.add(edge)
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
            if vMol.elementAt(i).isAvailable:
                num += 1
            i += 1
        return num

    def chooseVertex(self, vMol):
        num_atoms = len(vMol)
        index = -1
        min_num_linkages = 10
        i = 0
        while i < num_atoms:
            if vMol.elementAt(i).isAvailable:
                if vMol.elementAt(i).num_linkages < min_num_linkages:
                    min_num_linkages = vMol.elementAt(i).num_linkages
                    index = i
            i += 1
        return index

    def chooseSingleLinkedVertex(self, vMol):
        num_atoms = len(vMol)
        index = -1
        i = 0
        while i < num_atoms:
            if vMol.elementAt(i).isAvailable:
                if vMol.elementAt(i).num_linkages == 1:
                    return i
            i += 1
        return index

    def removeLinkage(self, vMol, index):
        i = 0
        while i < 5:
            if vMol.elementAt(index).linkage[i] != -1:
                while j < 5:
                    if vMol.elementAt(vMol.elementAt(index).linkage[i]).linkage[j] == index:
                        vMol.elementAt(vMol.elementAt(index).linkage[i]).linkage[j] = -1
                        vMol.elementAt(vMol.elementAt(index).linkage[i]).num_linkages -= 1
                        break
                    j += 1
            i += 1

    def replaceLinkage(self, vMol, vertex, index_i, index_j):
        i = 0
        while i < 5:
            if vMol.elementAt(index_i).linkage[i] == vertex:
                vMol.elementAt(index_i).linkage[i] = index_j
            if vMol.elementAt(index_j).linkage[i] == vertex:
                vMol.elementAt(index_j).linkage[i] = index_i
            i += 1

    def getNumLinkages(self, atom):
        num_linkages = 0
        i = 0
        while i < 5:
            if atom.linkage[i] != -1:
                num_linkages += 1
            i += 1
        return num_linkages


    public void remove( int vertex, Vector<Atom> vMol, Vector<Edge> vEdges, Vector<Edge> vRings ){


    def removeCommonVertex(self, path_i, path_j):
        array_i = path_i.split(" ")
        array_j = path_j.split(" ")
        new_path = ""
        vertex_i = int()
        vertex_j = int()
        isCommonVertex = False
        i = 0
        while len(array_i):
            vertex_i = Integer.parseInt(array_i[i])
            isCommonVertex = False
            while len(array_j):
                vertex_j = Integer.parseInt(array_j[j])
                if vertex_i == vertex_j:
                    isCommonVertex = True
                    break
                j += 1
            if not isCommonVertex:
                new_path += array_i[i] + " "
            i += 1
        return new_path.trim()

    def mergePaths(self, path_i, path_j):
        array_i = path_i.split(" ")
        array_j = path_j.split(" ")
        new_path = ""
        vertex_i = int()
        vertex_j = int()
        isCommonVertex = False
        i = 0
        while len(array_i):
            vertex_i = Integer.parseInt(array_i[i])
            isCommonVertex = False
            while len(array_j):
                vertex_j = Integer.parseInt(array_j[j])
                if vertex_i == vertex_j:
                    isCommonVertex = True
                    break
                j += 1
            if not isCommonVertex:
                new_path += array_i[i] + " "
            i += 1
        if 0 != len(length):
            return new_path.trim() + " " + path_j
        else:
            return path_j

    def removeLargeRings(self, vRings, num_atoms):
        array = []
        i = 0
        while i < len(vRings):
            array = vRings.elementAt(i).path.split(" ")
            if len(array):
                vRings.removeElementAt(i)
                i -= 1
            i += 1

    def clearEdges(self, vEdges, num_edges):
        array = []
        i = 0
        while i < len(vEdges):
            array = vEdges.elementAt(i).path.split(" ")
            if len(array):
                vEdges.removeElementAt(i)
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
            index_i = Integer.parseInt(array[i])
            while j < num_atoms:
                index_j = Integer.parseInt(array[j])
                if i != j:
                    dist = getDistance(vMol.elementAt(index_i).R, vMol.elementAt(index_j).R)
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
            if not isRing(vMol, vRings.elementAt(i)):
                vRings.removeElementAt(i)
                i -= 1
            i += 1

    public void removeRedundantRings( Vector<Edge> vRings ){
    public void keepSmallestRings( Vector<Edge> vRings ){

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
                edge_i = vRings.elementAt(i)
                path_1 = edge_i.path
                array_1 = path_1.split(" ")
                while j < len(vRings):
                    edge_j = vRings.elementAt(j)
                    path_2 = edge_j.path
                    array_2 = path_2.split(" ")
                    while len(array_1):
                        while len(array_2):
                            if Integer.parseInt(array_1[s1]) == Integer.parseInt(array_2[s2]):
                                mergeRings = True
                                vRings.elementAt(i).path = mergePaths(path_1, path_2)
                                vRings.removeElementAt(j)
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

