#!/usr/bin/python3
import re
from sys import argv,maxsize,path
path.insert(0,r'../')
from JorG.loadsave import POSCARloader
import numpy as np

BRUTAL='Ni'
CUTOFF=12.0
NUMBER=3

def read_moments(rawTxt):
    moments={}
    energy =0.0
    ISTOBEREAD=False
    i=0
    for line in rawTxt:
        if re.search('^\s+$',line):
            continue
        elif re.search('^\s*#',line):
            continue
        elif re.search('^-+$',line):
            continue
        elif re.search('^\s*tot',line):
            ISTOBEREAD = False
            continue
        elif "magnetization (x)" in line:
            ISTOBEREAD = True
            continue
        elif ISTOBEREAD:
            data = line.split()
            moments[int(data[0])] = float(data[4])
        elif "energy  without entropy" in line:
            energy = float(line.split()[-1])
    return moments,energy

    # 'comment','directions','cell',
    # 'cellSymmetry','cellVolume',
    # 'cellCenter','cellAtoms','atomNames'
def read_flip(name):
    with open(name+'/OUTCAR','r') as outcar:
        moments,energy = read_moments(outcar.readlines())
    loader   = POSCARloader(name+'/POSCAR',spam=False)
    loader.parse()
    directions = loader()['directions']
    cell       = loader()['cell']
    return moments,directions,cell,energy

from itertools import product
def get_crystal(moments,directions,cell,energy):
    crystal = []
    for mul in product(range(-1,1),repeat=3):
        boost = np.dot(mul,directions)
        for i,atom in enumerate(cell):
            crystal.append([i+1,atom[0],atom[1]+boost])
    
    return crystal

data = []
for name in argv[1:]:
    data.append((read_flip(name)))

recordA = data[0]
crystal = get_crystal(*recordA,)

distances= set([])
for atom in recordA[2]:
    d = np.around(np.linalg.norm(atom[1]),decimals=2)
    if d<CUTOFF and atom[0] in BRUTAL:
        distances.add(d)
distances = np.array([d for d in distances])
distances = np.sort(distances)
distances = distances[1:NUMBER+1]

systemOfEquations = []
dE = []
for i,recordB in enumerate(data[1:]):
    dE.append(recordA[3]-recordB[3])
    systemOfEquations.append(np.zeros(NUMBER))
    flipped = []
    for a,b in zip(recordA[0],recordB[0]):
        scalar = recordA[0][a] * recordB[0][b]
        if (abs(scalar) > 1e-5 and scalar < 0.0
             and recordA[2][a-1][0] in BRUTAL
             and recordB[2][b-1][0] in BRUTAL):
            flipped.append(a)

    for f,(j,atom) in product(flipped,enumerate(crystal)):
        if atom[0] not in flipped and atom[1] in 'Ni':
            d = np.around(np.linalg.norm(atom[2]-recordB[2][f-1][1]),decimals=2)
            if d in distances:
                systemOfEquations[i][np.argwhere(distances==d)] += recordA[0][atom[0]]*recordB[0][f]

systemOfEquations = np.array(systemOfEquations)
dE                = np.array(dE)

if systemOfEquations.shape[0] == systemOfEquations.shape[1]:
    sol = np.linalg.solve(systemOfEquations,dE)
else:
    sol = np.linalg.lstsq(systemOfEquations,dE,rcond=None)[0]

print("J in meV")
print(*(500*sol))
