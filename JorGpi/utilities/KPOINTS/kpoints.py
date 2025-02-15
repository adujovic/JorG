import numpy as np
from JorGpi.POSCARloader import POSCARloader

class KPOINTS:
    def __init__(self,POSCAR="POSCAR",resolution=0.001):
        self.multipliers = np.arange(1.0,89,resolution)
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

    @staticmethod
    def print_all(records,display):
        for record,error in records:
            print("{:3d} {:3d} {:3d} +/- {:4.3f}  {:7d}".format(*record,error,np.prod(record)).center(display))
