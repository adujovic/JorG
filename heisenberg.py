# -*- coding: utf-8 -*-

from sys import path
path.insert(0,r'../../')
import numpy as np
from itertools import product

class EquationSolver:
    def __init__(self,equations,vector):
        self.equations = np.array(equations)
        self.vector    = np.array(vector)
        self.solution  = None

    def solve(self,**kwargs):
        if self.solution is not None:
            return self.solution
        self.solve_system_of_equations(**kwargs)
        return self.solution

    def solve_system_of_equations(self, **kwargs):
        if len(self.equations) == 1:
            self.solution = self.vector[0]/self.equations[0]
        elif self.equations.shape[0] == self.equations.shape[1]:
            self.solution = np.linalg.solve(self.equations,self.vector,**kwargs)
        else:
            try:
                self.solution = np.linalg.lstsq(self.equations,self.vector,rcond=None,**kwargs)[0]
            except np.linalg.LinAlgError as err:
                print("%s: see:"%str(err))
                print(self.equations.shape,self.vector.shape)
                print(self.equations)
                print(self.vector)

    def remove_tautologies(self, **kwargs):
        # removing 0 = 0 equations !
        scale       = np.average(np.abs(self.equations))
        tautologies = np.argwhere(
                        np.apply_along_axis(
                          np.linalg.norm,1,self.equations)/scale<1e-3).flatten()
        self.equations = np.delete(self.equations,tuple(tautologies),axis=0)
        self.vector    = np.delete(self.vector   ,tuple(tautologies))
        return self.equations,self.vector

    @staticmethod
    def is_to_be_removed(i,j,equations):
        if i == j:
            return False
        scale = np.sum(np.abs(equations))
        inner_product = np.inner(equations[i],equations[j])
        norm_i = np.linalg.norm(equations[i])
        norm_j = np.linalg.norm(equations[j])
        if np.abs(inner_product - norm_j * norm_i)/scale < 1E-8:
            return True
        return False

    def remove_linears(self, **kwargs):
        # Based on https://stackoverflow.com/questions/
        #                  28816627/how-to-find-linearly-independent-rows-from-a-matrix
        # We remove lineary dependent rows
        remover=set([])
        for i,j in product(range(self.equations.shape[0]),repeat=2):
            if i in remover:
                continue
            if self.is_to_be_removed(i,j,self.equations):
                remover.add(j)
        if remover:
            self.equations = np.delete(self.equations,tuple(remover),axis=0)
            self.vector    = np.delete(self.vector   ,tuple(remover))
            return remover

    def remove_linear_combinations(self, secondSystem):
        scale            = np.sum(np.abs(self.equations))
        resultantSystem  = np.array([self.equations[0]])
        gram_schmidt     = np.array([self.equations[0]])
        self.equations   = np.delete(self.equations, (0), axis=0)
        resultantSystem2 = np.array([secondSystem[0]])
        secondSystem     = np.delete(secondSystem, (0), axis=0)
        while len(resultantSystem) < self.equations.shape[1]:
            if self.equations.size == 0:
                print("ERROR! Not enough equations! Please rerun.")
                exit(-3)
            tmpVector = np.copy(self.equations[0])
            for vector in gram_schmidt:
                tmpVector -= np.inner(tmpVector,vector)*vector/np.inner(vector,vector)
            if np.linalg.norm(tmpVector)/scale > 1e-4:
                gram_schmidt     = np.vstack((gram_schmidt,tmpVector))
                resultantSystem  = np.vstack((resultantSystem,self.equations[0]))
                self.equations   = np.delete(self.equations, (0), axis=0)
                resultantSystem2 = np.vstack((resultantSystem2,secondSystem[0]))
                secondSystem     = np.delete(secondSystem, (0), axis=0)
        self.equations = resultantSystem
        return resultantSystem,resultantSystem2

class DefaultMoments:
    def __init__(self,crystal):
        self.moments = [ atom[2] for atom in crystal ]
        self.moments.insert(0,0.0)
    def __call__(self,idx=0):
        return {'energy': 0.0, 'moments': self.moments}

class NaiveHeisenberg:
    def __init__(self,flippings,
                 crystal,crystal8):
        self.MAGMOMs   = DefaultMoments(crystal)
        self.flippings = flippings
        self.crystal   = crystal
        self.crystal8  = crystal8

    @staticmethod
    def name_interaction(first,second,distance):
        return "%s-%s @%.2f"%(first,second,distance)

    def initialize(self,mask,flipper):
        self.flipper           = np.array(flipper)
        self.mask              = mask
        self.numberOfElements  = self.mask.count('$')
        mul                    = self.numberOfElements*(self.numberOfElements + 1)//2
        self.systemOfEquations = np.zeros((len(self.flippings),len(flipper)*mul))
        self.interactionNames  = [ None ]*len(flipper)*mul
        self.avgMoments        = [  0.0 ]*len(flipper)*mul
        self.occMoments        = [  0   ]*len(flipper)*mul
        
        return mul

    def generate(self,mask,flipper):
        mul = self.initialize(mask,flipper)
        for (row,config),(I,atomI),atomJ in product(enumerate(self.flippings),
                                                    enumerate(self.crystal  ),
                                                              self.crystal8 ):
            if config[I] == config[atomJ[3]]:
                continue
            distance,offset = self.check_if_contributes(atomI,atomJ)
            if not distance:
                continue
            j = np.argwhere(np.abs(self.flipper - distance)<1e-2)
            if j.size: # geometric
                column = mul*j[0][0]+offset
                moment = np.abs(self.MAGMOMs()['moments'][I+1]\
                          *self.MAGMOMs(row+1)['moments'][atomJ[3]+1])
                self.systemOfEquations[row][column] += moment 
                self.avgMoments[column]             += moment
                self.occMoments[column]             += 1
                if self.interactionNames [column] is None:
                    self.interactionNames[column] = \
                        NaiveHeisenberg.name_interaction(atomI[0],atomJ[0],distance)
        self.clean_SoE()
        return self.systemOfEquations

    def clean_SoE(self):
        remover = [ i for i,name in enumerate(self.interactionNames) if name is None ]
        self.interactionNames  = [ name for name in self.interactionNames if name is not None ]
        self.interactionNames  = self.interactionNames[:len(self.flipper)]
        self.systemOfEquations = np.delete(self.systemOfEquations,remover,axis=1)
        self.systemOfEquations = self.systemOfEquations[:,0:len(self.flipper)]
        self.avgMoments        = np.delete(np.array(self.avgMoments),remover)
        self.occMoments        = np.delete(np.array(self.occMoments),remover)
        self.avgMoments        = np.array(self.avgMoments)/np.array(self.occMoments)
        self.avgMoments        = self.avgMoments[:len(self.flipper)]

    def offset_from_mask(self,symbol):
        element = self.mask.find(symbol+'$')
        return self.mask.count('$',0,element)

    def combine_offsets(self,offsetA,offsetB):
        return self.numberOfElements         *(self.numberOfElements         -1)//2 \
            - (self.numberOfElements-offsetA)*(self.numberOfElements-offsetA -1)//2 + offsetB


    def check_if_contributes(self,atomI,atomJ):
        if (atomI[0] not in self.mask
        or  atomJ[0] not in self.mask):
            return False,None
        offsetI = self.offset_from_mask(atomI[0])
        offsetJ = self.offset_from_mask(atomJ[0])
        if offsetI <= offsetJ:
            offset = self.combine_offsets(offsetI,offsetJ)
        else:
            offset = self.combine_offsets(offsetJ,offsetI)
        return np.round(np.linalg.norm(atomI[1]-atomJ[1]),2),offset

    def __str__(self):
        output = ""
        for I,i in product(range(len(self.MAGMOMs()['moments'])),
                           range(len(self.MAGMOMs))):
            if i == 0:
                output+="\n"
            output+="% .2f  "%self.MAGMOMs(i)['moments'][I+1]
        output+="\n"
        return output

def apply_mirrorsXYZ(dimensions,cell,reference=0):
    outputCell = []
    for p in product([-1,0,1],repeat=3):
        projection = np.array([p])
        translation = np.dot(projection,dimensions)[0]
        for i,atom in enumerate(cell):
            outputCell.append([atom[0],atom[1]+translation,atom[2],i])
    return outputCell
