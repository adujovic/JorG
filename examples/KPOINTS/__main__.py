from sys import path,argv
path.insert(0,r'../../')
import numpy as np
from POSCARloader import POSCARloader

class KPOINTS:
    def __init__(self,POSCAR="POSCAR",resolution=0.001):
        self.multipliers = np.arange(1.0,100,resolution)
        loader = POSCARloader(POSCAR)
        loader.parse()
        self.directions = loader()['directions']
        self.found      = []
        self.volume = np.linalg.det(self.directions)
        self.invert_directions()

    def invert_directions(self):
        self.invdirections = []
        for i in range(3):
            j = (i+1)%3
            k = (i+2)%3
            self.invdirections.append(np.cross(self.directions[j],self.directions[k])/self.volume)
        self.invdirections = np.array(self.invdirections)
        self.invdirections /= np.max([np.linalg.norm(b) for b in self.invdirections])

    def __call__(self,treshold):
        newlyFound = []
        for mul in self.multipliers:
            local = mul*self.invdirections
            err = np.sum(np.abs([np.linalg.norm(l)- np.int(np.linalg.norm(l)) for l in local]))
            if err >= treshold:
                continue
            justFound = [ np.int(np.linalg.norm(l)) for l in local ]
            if justFound not in self.found and 0 not in justFound:
                self.found.append(justFound)
                newlyFound.append([justFound,err])
        return newlyFound

def print_all(records,display):
    for record,error in records:
        print("{:3d} {:3d} {:3d} +/- {:4.3f}".format(*record,error).center(display))

def check(points,display):
    print("KPOINT multiplier proposition:".center(display))
    print("Good:".center(display))
    print_all(points(np.sqrt(1e-3)),display)
    print("Decent:".center(display))
    print_all(points(np.sqrt(1e-2)),display)
    print("Plausible:".center(display))
    print_all(points(np.sqrt(1e-1)),display)
    print("Borderline terrible:".center(display))
    print_all(points(np.sqrt(0.5)),display)
    print("Just don't:".center(display))
    print_all(points(np.sqrt(1.0)),display)

if __name__ == '__main__':
    display = 42
    for arg in argv[1:]:
        points = KPOINTS(arg)
        print(arg)
        check(points,display)
