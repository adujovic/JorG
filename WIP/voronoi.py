#!/usr/bin/python3

import numpy as np
import scipy.spatial #Voronoi, ConvexHull
from sys import path
path.insert(0,'../../')
from matplotlib import colors as mplColors
from itertools import product
from geometry.showcell import showCell

class Geometry:
    @staticmethod
    def radius_from_volume(volume):
        return np.power(0.75*volume/np.pi,np.reciprocal(3.0))

    @staticmethod
    def distance_to_tetragon(point,tetragon):
        dists = []
        vectors = np.array([ tetragon[i] - tetragon[0] for i in range(1,4) ])
        normal = np.zeros(4)
        for i in range(4):
            normal[i] = np.linalg.det(np.insert(vectors,0,np.eye(4)[i],axis=0))
        if np.linalg.norm(normal)<1e-5:
#            print(tetragon)
#            print(vectors)
#            print(np.insert(vectors,0,np.eye(4)[i],axis=0))
#            print(normal)
#            exit()
            return None
        else:
            normal /= np.linalg.norm(normal) 
        d = np.abs(np.dot(normal,point-tetragon[0]))
        if d < 0.1:
            return None
        return np.abs(np.dot(normal,point-tetragon[0]))

class Voronoi:
    multipliers = [0,-1,1]
    def __init__(self,data=None):
        self.data = data
        self.plotter = showCell(resolution=7)

        self.fill_if_read()

    def fill_if_read(self):
        if self.data is None:
            pass
        self.basis     = self.data['directions']
        self.basis     = np.append(self.basis,np.zeros((3,1)),axis=1)
        self.basis     = np.append(self.basis,[np.zeros(4)],axis=0)
        self.basis[-1,-1] = 1.0
        self.superCell = np.array([atom[1] for atom in self.data['cell']])
        self.atomNames = [atom[0] for atom in self.data['cell']]

        self.atomColors = {}
        self.atomWeight = {'Mn':0.4,'O': 0.6}
        for name in set(self.atomNames):
            self.atomColors[name] = mplColors.rgb2hex(0.8*np.random.rand(3))
    #        self.atomWeight[name] = (2*np.random.rand())**3
        self.center = np.mean(self.superCell, axis=0)
        self.center = np.append(self.center, 0.0)
        self.cutOff = np.linalg.norm(np.sum(self.basis,axis=0))

        self.points = [] #[np.append(atom,self.atomWeight[name]) for atom,name in zip(self.superCell,self.atomNames)]
        self.names  = [name for name in self.atomNames]

        self.data = None

    @staticmethod
    def support(center,radius,direction):
        project = np.dot(center,direction)
        return [project + radius, project - radius]

    @staticmethod
    def supportXYZ(center, radius):
        differences = []
        for v in np.identity(4):
            differences.append(Voronoi.support(center,radius,v))
        differences = np.array(differences)
        return [np.min(differences),np.max(differences)]

    def add_points(self,multipliers):
        for name,atom in zip(self.atomNames,self.superCell):
            vector = np.append(atom,self.atomWeight[name]) + np.dot(multipliers,self.basis)
#            vector = np.append(atom,self.atomWeight[name]) + np.dot(multipliers,self.basis)+ multipliers[-1]*self.atomWeight[name]*np.eye(4)[-1]
            self.add_point(name,vector)

    def add_point(self,name,vector):
        if(np.linalg.norm(vector-self.center) <= self.cutOff):
            self.points.append(vector)
            self.names.append(name)

    def fill_points(self):
        if self.data is None:
            pass
        else:
            self.fill_if_read()
#        for multipliers in product(self.multipliers,repeat=4):
        for multipliers in product(self.multipliers,self.multipliers,self.multipliers,[0,-2,-1,10]):
            self.add_points(multipliers)
        self.points = np.array(self.points)

    def get_Voronoi_diagram(self, name=None,
                                  save=False):
        self.fill_points()
        self.diagram = scipy.spatial.Voronoi(self.points)
        self.regions = [self.diagram.regions[
                        self.diagram.point_region[i]]
                        for i,atom in enumerate(self.superCell)]
        self.calculate_all()
        np.asarray(self.WignerSeitzRadia)
        if name is not None:
            np.savetxt(name, self.WignerSeitzRadia)
        elif save:
            np.savetxt("WignerSeitzRadia.dat", self.WignerSeitzRadia)

    def clear_vertices(self):
        for i,region in enumerate(self.regions):
            infty = np.where(np.array(region)<0)
            self.regions[i] = np.delete(region,*infty)

    def add_convex_hull(self,vertices,atom,name):
        for simplex in self.convexHull.simplices:
            distance = Geometry.distance_to_tetragon(atom,vertices[simplex])
            if distance is None:
                continue
            if distance < self.radius:
                self.radius = distance
                self.plotter.add_polygon(vertices[simplex][:,:-1],
                                     alpha=0.4,
                                     color=self.atomColors[name])

    def calculate_all(self):
        self.clear_vertices()
        self.aspect_data      = []
        self.WignerSeitzRadia = []
        for name,atom,region in zip(self.atomNames,
                                    self.superCell,
                                    self.regions):
            self.WignerSeitzRadia.append(
                    self.calculate_radia(name,np.append(atom,self.atomWeight[name]),region))
        self.plotter.set_aspect(self.aspect_data)

    def calculate_radia(self,name,atom,region):
        vertices = np.array([self.diagram.vertices[indx] for indx in region])
        self.radius = self.cutOff*2*10
        try:
            self.convexHull = scipy.spatial.ConvexHull(vertices)
        except scipy.spatial.qhull.QhullError:
            print(vertices)
            exit()
        WSradius = Geometry.radius_from_volume(self.convexHull.volume)
        self.add_convex_hull(vertices,atom,name)
        self.plotter.add_sphere(atom,self.radius,color=self.atomColors[name],alpha=1.0)
        self.aspect_data.append(Voronoi.supportXYZ(atom,self.radius))
        return (self.radius,WSradius)

    def show(self,name=None):
        self.plotter.show(name)
