#!/usr/bin/python3

import numpy as np
import scipy.spatial #Voronoi, ConvexHull
from sys import argv,path
path.insert(0,'../')
from matplotlib import colors as mplColors
from itertools import product
from voronoi.showcell import showCell

class Voronoi:
    multipliers = [-1,0,1]
    def __init__(self,data=None):
        self.data = data
        self.plotter = showCell(resolution=7)

        self.fill_if_read()

    def fill_if_read(self):
        if self.data is None:
            pass
        self.basis = self.data['directions']
        self.superCell = np.array([atom[1] for atom in self.data['cell']])
        self.atomNames = [atom[0] for atom in self.data['cell']]
    
        self.atomColors = {}
        for name in set(self.atomNames):
            self.atomColors[name] = mplColors.rgb2hex(0.8*np.random.rand(3))
        self.center = np.mean(self.superCell, axis=0)
        self.cutOff = np.linalg.norm(np.sum(self.basis,axis=0))

        self.points = []
        self.names  = []
        for name,atom in zip(self.atomNames,self.superCell):
            self.points.append(atom)
            self.names.append(name)
    
        self.data = None

    @staticmethod
    def support(center,radius,direction):
        project = np.dot(center,direction)
        return [project + radius, project - radius]
    
    @staticmethod
    def supportXYZ(center, radius):
        differences = []
        for v in np.identity(3):
            differences.append(Voronoi.support(center,radius,v))
        differences = np.array(differences)
        return [np.min(differences),np.max(differences)]

    def fill_points(self):
        if self.data is None:
            pass
        else:
            self.fill_if_read()
        for i,j,k in product(self.multipliers,repeat=3):
            if i == 0 and j == 0 and k == 0:
                continue
            for name,a in zip(self.atomNames,self.superCell):
                v = a + i*self.basis[0]+j*self.basis[1]+k*self.basis[2]
                if(np.linalg.norm(v-self.center) <= self.cutOff):
                    self.points.append(v)
                    self.names.append(name)
        self.points = np.array(self.points)

    def get_Voronoi_diagram(self):
        self.diagram = scipy.spatial.Voronoi(self.points)
        self.regions = [self.diagram.regions[
                         self.diagram.point_region[i]]
                        for i,atom in enumerate(self.superCell)]

    def clear_vertices(self):
        for i,region in enumerate(self.regions):
            infty = np.where(np.array(region)<0)
            self.regions[i] = np.delete(region,*infty)

    def show(self):
        self.clear_vertices()
        aspect_data = []
        for name,atom,region in zip(self.atomNames,
                                    self.superCell,
                                    self.regions):
            vertices = np.array([self.diagram.vertices[indx] for indx in region])
            radius = self.cutOff*2;
            convexHull = scipy.spatial.ConvexHull(vertices)
            volume = convexHull.volume
            WSradius = np.power(3*volume/(4*np.pi),1.0/3.0)
            for simplex in convexHull.simplices:
                normal = np.cross(vertices[simplex[1]]-vertices[simplex[0]],vertices[simplex[2]]-vertices[simplex[0]])
                normal /= np.linalg.norm(normal)
                distance = np.abs(np.dot(normal,atom-vertices[simplex[0]]))
                if distance < radius:
                    radius = distance
                self.plotter.add_polygon([region[simplex]], alpha=0.2)
            self.plotter.add_sphere(atom,radius,color=self.atomColors[name],alpha=1.0)
            aspect_data.append(Voronoi.supportXYZ(atom,radius))
            print(name,radius,WSradius)
        self.plotter.set_aspect(aspect_data)
        self.plotter.show()
