#!/usr/bin/env python
import sys
sys.dont_write_bytecode = True

class Fragment(object):
    linked_atoms = str()
    type_ = int()

    #  1: ring, 2: substituent, 3: linker, 4: functional group, 5: non-functional group
    def __init__(self):
        self.linked_atoms = ""
        self.type_ = 0

