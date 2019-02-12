#!/usr/bin/env python
import sys
sys.dont_write_bytecode = True

class LinkedAtomGroup(object):

    def __init__(self):
        self.linked_atoms = ""
        self.linked_bonds = ""
        self.i = self.j = -1