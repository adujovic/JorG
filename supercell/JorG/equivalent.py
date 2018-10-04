import numpy as np
from aux.periodic import maskFull

def find_unique_flips(referenceAtom,crystal,symmetry,cutOff,atomTypeMask=maskFull,logAccuracy=2):
    distances = [ ]
    for atom in crystal:
       d = np.around(np.linalg.norm(atom[1]-referenceAtom[1]),decimals=logAccuracy )
       if (d,atom[0]) not in distances and d > 0.0 and d <= cutOff:
           distances.append((d,atom[0]))
    
    distances = np.array(distances, dtype=[('distance', np.float), ('element', 'U3')])
    distances.sort(order='distance')
    
                                        #for d,n in distances:
                                        #    print(n,d)
                                        #for atom in crystal:
                                        #  print("%s %f %f %f"%(atom[0],atom[1][0],atom[1][1],atom[1][2]))
    
    checker = [False for x in range(len(crystal))] 
    flipper = [[False,[]] for x in range(len(crystal))] 
    case = 1
    for (d,n) in distances:
        if n in atomTypeMask:
            for i,atom in enumerate(crystal):
                if not checker[i]:
                    if d == np.around(np.linalg.norm(atom[1]-referenceAtom[1]),decimals=logAccuracy) and atom[0] == n:
                        for j in np.argwhere(symmetry['equivalent_atoms'] == symmetry['equivalent_atoms'][i]).flatten():
                            if d == np.around(np.linalg.norm(crystal[j][1]-referenceAtom[1]),decimals=logAccuracy):
                                print("Case no: %d, for %s with d= %f"%(case,crystal[j][0],d))
                                checker[j] = True
                        if checker[i]:
                            flipper[i][0] = True
                            flipper[i][1] = atom
                            case += 1
    return flipper

