#!/usr/bin/env python
import sys
sys.dont_write_bytecode = True

class Atom(object):
    atom_name = str()
    element = str()
    line = str()
    R = []
    isAvailable = bool()
    isRingAtom = bool()
    isAromatic = bool()
    isBridgingAtom = bool()
    ringIndex = []
    numRingAtoms = []

    #  e.g., five-membered ring, six-membered ring
    ringHasCabonyl = bool()
    ringHasAmide = bool()
    isCarbonyl = bool()
    isPyrrole = bool()
    isNeutralPyridine = bool()
    isProtonatedPyridine = bool()
    isProtonatedPyrimidine = bool()
    isFuran = bool()
    isFuranose = bool()
    isProtonatedNitrogen = bool()
    isDeprotonatedSulfur = bool()
    isDeprotonatedOxygen = bool()
    isImineCarbon = bool()
    isProtonatedImineGroup = bool()
    isKetone = bool()
    isPOxide = bool()
    isSOxide = bool()
    isAmide = bool()
    isEster = bool()
    isThiocarbonylS = bool()
    numHydrogenAtoms = int()
    numNitrogenAtoms = int()
    numCarbonAtoms = int()
    numOxygenAtoms = int()
    numSulfurAtoms = int()
    isMethyl = bool()
    isConjugated = bool()
    linkage = []
    bondType = []
    bondOrder = []
    totalBondOrder = int()

    #  Bond Types:
    #  1 = single
    #  2 = double
    #  3 = triple
    #  am = amide
    #  ar = aromatic
    #  du = dummy
    #  un = unknown (cannot be determined from the parameter tables)
    #  nc = not connected
    num_linkages = int()
    linked_atoms = str()
    atomType = str()
    atomPriority = int()
    charge = float()
    group = int()
    formal_charge = float()
    charge_calc = float()
    charge_temp = float()

    def __init__(self):
        self.atom_name = None
        self.element = None
        self.line = None
        self.R = [None]*3
        self.isAvailable = True
        self.isRingAtom = False
        self.isAromatic = False
        self.isBridgingAtom = False
        self.ringIndex = [None]*3
        self.numRingAtoms = [None]*3
        i = 0
        while i < 3:
            self.ringIndex[i] = -1
            self.numRingAtoms[i] = 0
            i += 1
        self.ringHasCabonyl = False
        self.ringHasAmide = False
        self.isCarbonyl = False
        self.isPyrrole = False
        self.isNeutralPyridine = False
        self.isProtonatedPyridine = False
        self.isProtonatedPyrimidine = False
        self.isFuran = False
        self.isFuranose = False
        self.isProtonatedNitrogen = False
        self.isDeprotonatedSulfur = False
        self.isDeprotonatedOxygen = False
        self.isImineCarbon = False
        self.isProtonatedImineGroup = False
        self.isKetone = False
        self.isPOxide = False
        self.isSOxide = False
        self.isAmide = False
        self.isEster = False
        self.isThiocarbonylS = False
        self.numHydrogenAtoms = 0
        self.numNitrogenAtoms = 0
        self.numCarbonAtoms = 0
        self.numOxygenAtoms = 0
        self.numSulfurAtoms = 0
        self.isMethyl = False
        self.isConjugated = False
        self.linkage = [None]*5
        self.bondType = [None]*5
        self.bondOrder = [None]*5
        i = 0
        while i < 5:
            self.linkage[i] = -1
            self.bondType[i] = None
            self.bondOrder[i] = 0
            i += 1
        self.totalBondOrder = 0
        self.num_linkages = 0
        self.linked_atoms = None
        self.atomType = None
        self.atomPriority = -1
        self.charge = 0
        self.group = -1
        self.formal_charge = 0
        self.charge_calc = 0
        self.charge_temp = 0

