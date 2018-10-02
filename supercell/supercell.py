#!/usr/bin/python3

import numpy as np
import spglib as sp
import re
from itertools import permutations as per

cutOff = 20.0

comment = ""
directions = []
normals = np.identity(3)
cell = []
cellVolume = 0.0
cellCenter = np.zeros(3)
cellAtoms = np.array([],dtype=np.int)
atomNames = []
cellSize = 1
cellInputType = 'x'
vectors = False

""" Templates:
    according to vasp.wiki:
          direct coords are for 6th lines starting with 'D' or 'd'
      carthesian coords are for 6th lines starting with 'C', 'c', 'K' or 'k' """

directTemplate     = "Dd"
carthTemplate      = "CcKk"
selectiveTemplate  = "Ss"
                   
""" Reading POSCAR file.
      TODO: bulletproofing """

offset = 0

with open("POSCAR","r+") as inFile:
    for i,raw in enumerate(inFile.readlines()):
        line = re.sub("^\s*","",raw)  # remove all blank characters from begining of the line
        line = re.sub("\s+"," ",line) # replace all blank characters in a row to a single space
        if i == 0:
            comment = line[:-1]
        elif i == 1:
            scale = np.float64(line[:-1]) # scaling factor
            if(scale < 0.0):
                scale = 1.0
        elif i in range(2,5):
            directions.append(scale*np.fromstring(line,sep=" ")) # crystal directions
            cellCenter += 0.5*directions[-1]
        elif i == 5:
            cellVolume = np.abs(np.linalg.det(np.array(directions)))
            if re.match("\D",line):
                offset += 1
                atomNames=line[:-1].split(" ")
            else:    
                cellAtoms = np.fromstring(line,sep=" ",dtype=np.int)
                cellSize = np.sum(cellAtoms)
        if i == 6 and offset == 1:
            cellAtoms = np.fromstring(line,sep=" ",dtype=np.int)
            cellSize = np.sum(cellAtoms)
        elif i == 6 + offset:
            cellInputType = line[0]
            if cellInputType in carthTemplate:
                vectors = normals
            elif cellInputType in directTemplate:
                vectors = directions
            elif cellInputType in selectiveTemplate:
                offset += 1
            else:
                print("Error in POSCAR: unknown input in line %d: %s"%(i,line[:-1]))
                exit(-1)
        elif not vectors and i == 7:
            cellInputType = line[0]
            if cellInputType in carthTemplate:
                vectors = normals
            elif cellInputType in directTemplate:
                vectors = directions
            else:
                print("Error in POSCAR: unknown input in line %d: %s"%(i,line[:-1]))
                exit(-1)
        elif i in range(7+offset,7+offset+cellSize):
            found = re.search("([\-\+]?\d+\.?\d*)\s([\-\+]?\d+\.?\d*)\s([\-\+]?\d+\.?\d*)",line)
            if found:
                atomType = np.flatnonzero(cellAtoms)
                cellAtoms[atomType[0]] -= 1
                atom  = [atomType[0],np.zeros(3)]
                coords = np.fromstring(found.group(0),sep=" ")
                for coord,vector in zip(coords,vectors):
                    atom[1] += coord*vector 
                cell.append(atom)    
            found = re.search("[\-\+]?\d+\.?\d*\s[\-\+]?\d+\.?\d*\s[\-\+]?\d+\.?\d*\s([a-zA-Z]+)",line)
            if found:
                if found.group(1) not in atomNames:
                    atomNames.append(found.group(1))

multiplyers = []

bestScore = 0.0
bestPermutation = (0,1,2) 
for p in per(range(3)):
    score = 1.0
    for i,d in enumerate(directions):
        score *= d[p[i]]
    if np.abs(score) > bestScore:
        bestScore = np.abs(score)
        bestPermutation = p

for d,p in zip(directions,
               bestPermutation):
    multiplyers.append(1+int(cutOff/d[p]))

periodicTable = {0  :"H", 1  :"He", 2  :"Li", 3  :"Be", 4  :"B", 5  :"C", 6  :"N", 7  :"O", 8  :"F", 9 :"Ne", 10 :"Na", 11 :"Mg", 12 :"Al", 13 :"Si", 14 :"P", 15 :"S", 16 :"Cl", 17 :"Ar", 18 :"K", 19 :"Ca", 20 :"Sc", 21 :"Ti", 22 :"V", 23 :"Cr", 24 :"Mn", 25 :"Fe", 26 :"Co", 27 :"Ni", 28 :"Cu", 29 :"Zn", 30 :"Ga", 31 :"Ge", 32 :"As", 33 :"Se", 34 :"Br", 35 :"Kr", 36 :"Rb", 37 :"Sr", 38 :"Y", 39 :"Zr", 40 :"Nb", 41 :"Mo", 42 :"Tc", 43 :"Ru", 44 :"Rh", 45 :"Pd", 46 :"Ag", 47 :"Cd", 48 :"In", 49 :"Sn", 50 :"Sb", 51 :"Te", 52 :"I", 53 :"Xe", 54 :"Cs", 55 :"Ba", 56 :"La", 57 :"Ce", 58 :"Pr", 59 :"Nd", 60 :"Pm", 61 :"Sm", 62 :"Eu", 63 :"Gd", 64 :"Tb", 65 :"Dy", 66 :"Ho", 67 :"Er", 68 :"Tm", 69 :"Yb", 70 :"Lu", 71 :"Hf", 72 :"Ta", 73 :"W", 74 :"Re", 75 :"Os", 76 :"Ir", 77 :"Pt", 78 :"Au", 79 :"Hg", 80 :"Tl", 81 :"Pb", 82 :"Bi", 83 :"Po", 84 :"At", 85 :"Rn", 86 :"Fr", 87 :"Ra", 88 :"Ac", 89 :"Th", 90 :"Pa", 91 :"U", 92 :"Np", 93 :"Pu", 94 :"Am", 95 :"Cm", 96 :"Bk", 97 :"Cf", 98 :"Es", 99:"Fm", 100:"Md", 101:"No", 102:"Lr", 103:"Rf", 104:"Db", 105:"Sg", 106:"Bh", 107:"Hs", 108:"Mt", 109:"Ds", 110:"Rg", 111:"Cn", 112:"Nh", 113:"Fl", 114:"Mc", 115:"Lv", 116:"Ts", 117:"Og"}
crystal = []

numberOfAtoms = 0
output = ""
for x in range(-multiplyers[0],multiplyers[0]+1):
    for y in range(-multiplyers[1],multiplyers[1]+1):
        for z in range(-multiplyers[2],multiplyers[2]+1):
            for atom in cell:
                position = np.copy(atom[1])
                for a,n in zip([x,y,z],directions):
                    position += a*n
                if np.linalg.norm(position - cellCenter) < cutOff :    
                    numberOfAtoms += 1
                    if atom[0] < len(atomNames):
                        output += "%s"%atomNames[int(atom[0])]
                    else:    
                        output += "%s"%periodicTable[int(atom[0])]
                    for p in position:
                        output += " %.10f"%p
                    output += "\n"

print(numberOfAtoms)
print("")
print(output)

exit()
print(comment)
print(directions)
print(cell)
