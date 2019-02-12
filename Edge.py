#!/usr/bin/env python
import sys
sys.dont_write_bytecode = True

class Edge(object):
    i = int()
    j = int()
    path = str()

    def __init__(self):
    	self.path = ""
        self.i = self.j = -1

