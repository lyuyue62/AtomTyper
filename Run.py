#!/usr/bin/python

import sys,os
import shlex, subprocess

errlist = open("error.txt","w")


for mol2 in os.listdir("mol2_modified"):
    cmd = "python2 AtomTyper.py mol2_modified/%s" % mol2
    proc = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if len(stderr) > 0:
        errlist.write("%s\n" % mol2)

errlist.close()

