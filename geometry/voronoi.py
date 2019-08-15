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
    def distance_to_triangle(point,triangle):
        normal = np.cross(triangle[1]-triangle[0],
                          triangle[2]-triangle[0])
        normal /= np.linalg.norm(normal)
        return np.abs(np.dot(normal,point-triangle[0]))

class Voronoi:
    multipliers = [-1,0,1]
    def __init__(self,data=None):
        self.data = data
        self.plotter = showCell(resolution=37)

        self.fill_if_read()

    def fill_if_read(self):
        if self.data is None:
            pass
        self.basis     = self.data['directions']
        self.superCell = np.array([atom[1] for atom in self.data['cell']])
        self.atomNames = [atom[0] for atom in self.data['cell']]

        self.atomColors = {}
        for name in set(self.atomNames):
            self.atomColors[name] = mplColors.rgb2hex(0.8*np.random.rand(3))
        self.center = np.mean(self.superCell, axis=0)
        self.cutOff = np.linalg.norm(np.sum(self.basis,axis=0))

        self.points = [atom for atom in self.superCell]
        self.names  = [name for name in self.atomNames]

        self.data = None

    @staticmethod
    def support(center,radius,direction):
        project = np.dot(center,direction)
        return [project + radius, project - radius]

    @staticmethod
    def support_xyz(center, radius):
        differences = []
        for versor in np.identity(3):
            differences.append(Voronoi.support(center,radius,versor))
        differences = np.array(differences)
        return [np.min(differences),np.max(differences)]

    def add_points(self,multipliers):
        for name,atom in zip(self.atomNames,self.superCell):
            vector = atom + np.dot(multipliers,self.basis)
            self.add_point(name,vector)

    def add_point(self,name,vector):
        if np.linalg.norm(vector-self.center) <= self.cutOff :
            self.points.append(vector)
            self.names.append(name)

    def fill_points(self):
        if self.data is None:
            pass
        else:
            self.fill_if_read()
        for multipliers in product(self.multipliers,repeat=3):
            self.add_points(multipliers)
        self.points = np.array(self.points)

    def get_voronoi_diagram(self, name=None,
                                  save=False):
        self.fill_points()
        self.diagram = scipy.spatial.Voronoi(self.points)
        self.regions = [self.diagram.regions[
                        self.diagram.point_region[i]]
                        for i,atom in enumerate(self.superCell)]
        self.calculate_all()
        np.asarray(self.wignerSeitzRadia)
        if name is not None:
            np.savetxt(name, self.wignerSeitzRadia)
        elif save:
            np.savetxt("wignerSeitzRadia.dat", self.wignerSeitzRadia)

    def clear_vertices(self):
        for i,region in enumerate(self.regions):
            infty = np.where(np.array(region)<0)
            self.regions[i] = np.delete(region,*infty)

    def add_convex_hull(self,vertices,atom,name):
        for simplex in self.convexHull.simplices:
            distance = Geometry.distance_to_triangle(atom,vertices[simplex])
            if distance < self.radius:
                self.radius = distance
            self.plotter.add_polygon(vertices[simplex],
                                     alpha=0.4,
                                     color=self.atomColors[name])

    def calculate_all(self):
        self.clear_vertices()
        self.aspect_data      = []
        self.wignerSeitzRadia = []
        for name,atom,region in zip(self.atomNames,
                                    self.superCell,
                                    self.regions):
            self.wignerSeitzRadia.append(
                    self.calculate_radia(name,atom,region))
        self.plotter.set_aspect(self.aspect_data)

    def calculate_radia(self,name,atom,region):
        vertices = np.array([self.diagram.vertices[indx] for indx in region])
        self.radius = self.cutOff*2
        self.convexHull = scipy.spatial.ConvexHull(vertices)
        WSradius = Geometry.radius_from_volume(self.convexHull.volume)
        self.add_convex_hull(vertices,atom,name)
        self.plotter.add_sphere(atom,self.radius,color=self.atomColors[name],alpha=1.0)
        self.aspect_data.append(Voronoi.support_xyz(atom,self.radius))
        return (self.radius,WSradius)

    def show(self,name=None):
        self.plotter.show(name)
