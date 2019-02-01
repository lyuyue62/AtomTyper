#!/usr/bin/env python

class Bond(object):
    i = int()
    j = int()
    order = int()
    bondType = str()

    #  *** For atomic partial charge assignment ***
    increment = float()
    index = int()

    def __init__(self):
        self.i = -1
        self.j = -1
        self.order = 0
        self.bondType = None
        self.increment = 0
        self.index = -1

