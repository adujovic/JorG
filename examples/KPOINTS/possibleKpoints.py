#!/usr/bin/python3
# -*- coding: utf-8 -*-

from sys import path
path.insert(0,r'../../')
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

def check(treshold,multipliers,found):
    for mul in multipliers:
        loc = mul*invdirections
        err = 0.0
        for l in loc:
            err += np.abs(np.linalg.norm(l)- np.int(np.linalg.norm(l)))
        if err < treshold:
            newlyFound = []
            for l in loc:
                newlyFound.append(np.int(np.linalg.norm(l)))
            if newlyFound not in found and 0 not in newlyFound:
                print("{:3d} {:3d} {:3d} +/- {:4.3f}".format(*newlyFound,err).center(display))
                found.append(newlyFound)
    return found

print("KPOINT multiplier proposition:".center(display))
found = []
multipliers = np.arange(1.0,100,resolution)
print("Good:".center(display))
found = check(np.sqrt(1e-3),multipliers,found)
print("Decent:".center(display))
found = check(0.1,multipliers,found)
print("Plausible:".center(display))
found = check(np.sqrt(1e-1),multipliers,found)
print("Borderline terrible:".center(display))
found = check(0.5,multipliers,found)
print("Just don't:".center(display))
found = check(1.0,multipliers,found)
