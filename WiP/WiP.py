# -*- coding: utf-8 -*-

from sys import argv,path
path.insert(0,r'../')
import numpy as np

#from subprocess import call
#import errno
#import re
#import spglib
#import shutil
#from itertools import product
#
#import JorG.symmetry
#from JorG.argv import options
#from JorG.format import color, print_case, print_vector
#from JorG.format import print_crystal, print_moments, print_label
#import JorG.loadsave as loadsave
#import JorG.generator as generator

class errors:
    failed_to_generate = -201
    systemerror = -1

from datetime import datetime
from os import makedirs
import JorG.loadsave as loadsave
class StreamHandler:
    def __init__(self,*args,**kwargs):
        self.streams = []
        for arg in args:
            name = StreamHandler.fix_name(arg)
            StreamHandler.make_directory(name)
            self.streams.append(name)

    @staticmethod
    def fix_name(name):
       if name is None:
           # if output directory is not given:
           name = "output/"+datetime.now().strftime("%Y%m%d%H%M%S")
       else:
           # remove multiple '/' and possible '/' at the end
           name = re.sub('/+','/',name)
           name = re.sub('/$','',name)
       return name

    @staticmethod
    def make_directory(name):
        # creating output path
        try:
            makedirs(name)
        except OSError as err:
            if err.errno == errno.EEXIST:
                print("%s exists: Data will be overwritten!"%name)

    def __call__(self, idx=0):
        return self.streams[idx]

    def load_VASP(self,POSCARfile,INCARfile):
        load_POSCAR          = loadsave.POSCARloader(POSCARfile)
        load_POSCAR.parse()
        readData             = load_POSCAR(0)
        oldMoments,incarData = loadsave.load_INCAR (readData['cell'],INCARfile,atomNames=readData['atomNames'])
        return readData,oldMoments,incarData

from copy import copy
class standard:
    values = {'scale'     : 5.0,
              'vector'    : 0.1,
              'background': [253,246,227],
              'radius'    : 10.0,
              'center'    : [0.0, 0.0, 0.0],
              'color'     : [176,255,176]}
    @staticmethod
    def fix(**kwargs):
        result = copy(standard.values)
        result.update(kwargs)
        return result

class JmolVisualization:
    scriptText =\
"""vibration ON
vector ON 
vector SCALE %.2f
vector %.2f
set bondRadiusMilliAngstroms 0 
background [%d,%d,%d]
isoSurface SURF center {%.5f %.5f %.5f} sphere %.5f
color $SURF translucent [%d,%d,%d]"""
    def create_script(outDirName,**kwargs):
        kwargs = standard.fix(**kwargs)
        with open(outDirName+"/script.jmol","w+") as stream:
            stream.write(JmolVisualization.scriptText%(
                kwargs['scale'],
                kwargs['vector'],
               *kwargs['background'],
               *kwargs['center'],
                kwargs['radius'],
               *kwargs['color'],))


from os import remove
class TemporaryFiles:
    names     = [".input",".supercell",".directions"]
    extension = ".dat"
    def __init__(self, suffix=str(np.random.randint(10000000,99999999)), prefix=''):
        self.suffix = suffix
        self.prefix = prefix
        for name in self.names:
            with open(self.prefix+name+self.suffix+self.extension,"w+") as w:
                w.write('hello')

    def write_input(self,allFlippable,crystal):
        name = self.prefix+self.names[0]+self.suffix+self.extension
        with open(name,"w+") as isingFile:
            for i in allFlippable:
                isingFile.write("%d %.8f %.8f %.8f %.2f\n"%(
                i,*crystal[i][1],crystal[i][2]))

    def write_supercell(self,crystal):
        name = self.prefix+self.names[1]+self.suffix+self.extension
        with open(name,"w+") as supercellFile:
            for i,atom in enumerate(crystal):
                supercellFile.write("%d %.8f %.8f %.8f %.2f\n"%(
                i,*atom[1],atom[2]))

    def write_directions(self,directions):
        name = self.prefix+self.names[2]+self.suffix+self.extension
        with open(name,"w+") as dirFile:
            for d in directions:
                dirFile.write("%.8f %.8f %.8f\n"%tuple(d))

    def __str__(self):
        iName = self.prefix+self.names[0]+self.suffix+self.extension
        sName = self.prefix+self.names[1]+self.suffix+self.extension
        dName = self.prefix+self.names[2]+self.suffix+self.extension
        return "%s %s %s"%(dName,sName,iName)

    def __del__(self):
        for name in self.names:
            try:
                remove(self.prefix+name+self.suffix+self.extension)
            except OSError:
                print("No file %s"%(self.prefix+name+self.suffix+self.extension))

from JorG.PeriodicTable import periodicTableElement
def from_refined(refinedCell):
    directions = np.array(refinedCell[0])
    cell = []
    for refinedAtom,atomType in zip(refinedCell[1],refinedCell[2]):
        newPosition = np.dot(refinedAtom,refinedCell[0])
        cell.append((periodicTableElement[atomType-1],newPosition))
    return cell,directions
