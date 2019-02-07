#!/usr/bin/env python
import sys
import math
import numpy as np


class SmallMolecule:

	def __init__(self):
		self.name = "LIG"
		self.charge = 0
		self.atoms = []
		self.bonds = []
		self.angles = []
		self.dihedrals = []
		self.impropers = []

	def setLinkage(self):
		for i in range(len(self.atoms)):
			self.atoms[i].num_linkages = 0
		
			for j in range(5):
				self.atoms[i].linkage[j] = -1
		
		
		for i in range(len(self.bonds)):
			for j in range(5):
				if  self.atoms[self.bonds[i].i].linkage[j] == -1 :
					self.atoms[self.bonds[i].i].linkage[j] = self.bonds[i].j
					self.atoms[self.bonds[i].i].bondType[j] = self.bonds[i].bondType
					self.atoms[self.bonds[i].i].num_linkages += 1			
					break
		
			for j in range(5):
				if  self.atoms[self.bonds[i].j].linkage[j] == -1:
					self.atoms[self.bonds[i].j].linkage[j] = self.bonds[i].i
					self.atoms[self.bonds[i].j].bondType[j] = self.bonds[i].bondType
					self.atoms[self.bonds[i].j].num_linkages += 1
					break
			
			if self.bonds[i].bondType.lower() ==  "ar":
				self.atoms[self.bonds[i].i].isAromatic = True
				self.atoms[self.bonds[i].j].isAromatic = True


	def setBridgingAtoms(self):
		num_rings = 0
		
		for i in range(len(self.atoms)):
			num_rings = 0
		
			for j in range(3):
				if self.atoms[i].ringIndex[j] != -1:
					num_rings += 1
		
			if num_rings >= 2 :
				self.atoms[i].isBridgingAtom = True
	