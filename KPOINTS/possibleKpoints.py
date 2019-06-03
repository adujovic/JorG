#!/usr/bin/python3
# -*- coding: utf-8 -*-

from sys import maxsize,path
path.insert(0,r'../')
import numpy as np
from JorG.loadsave import POSCARloader

display    = 42
resolution = 0.001

loader = POSCARloader("POSCAR")
loader.parse()
directions = loader()['directions'] 
volume = np.linalg.det(directions)

invdirections = []
for i in range(3):
    j = (i+1)%3
    k = (i+2)%3
    invdirections.append(np.cross(directions[j],directions[k])/volume)
invdirections = np.array(invdirections)
invdirections /= np.max([np.linalg.norm(b) for b in invdirections])

def check(treshold):
    global found,multipliers
    for mul in multipliers:
        loc = mul*invdirections
        err = 0.0
        for l in loc:
            test = np.linalg.norm(l)
            err += np.abs(test - np.int(test))
        if err < treshold:
            newlyFound = []
            for l in loc:
                newlyFound.append(np.int(np.linalg.norm(l)))
            if 0 in newlyFound:
                continue
            if newlyFound not in found:
                print("{:3d} {:3d} {:3d} +/- {:4.3f}".format(*newlyFound,err).center(display))
                found.append(newlyFound)

print("KPOINT multiplier proposition:".center(display))
found = []
multipliers = np.arange(1.0,100,resolution)
print("Good:".center(display))
check(np.sqrt(1e-3))
print("Decent:".center(display))
check(0.1)
print("Plausible:".center(display))
check(np.sqrt(1e-1))
print("Borderline terrible:".center(display))
check(0.5)
print("Just don't:".center(display))
check(1.0)

