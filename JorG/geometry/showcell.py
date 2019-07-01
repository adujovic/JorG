#!/usr/bin/python3

from matplotlib import colors as mplColors
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
import numpy as np
from sys import path
path.insert(0,'../../')


class showCell:
    figure  = None # matplotlib figure
    subplot = None # subplot

    def make_spherical(self):
        self.phi   = np.linspace(0, 2*np.pi, self.resolution) # spherical coordinates
        self.theta = np.linspace(0,   np.pi, self.resolution)

    def make_points(self,position,radius):
        return np.array([
                  position[0]+radius*np.outer(np.cos(self.phi), np.sin(self.theta)),
                  position[1]+radius*np.outer(np.sin(self.phi), np.sin(self.theta)),
                  position[2]+radius*np.outer(np.ones(np.size(self.phi)), np.cos(self.theta))])
    
    def __init__(self,resolution = 13, alpha = 0.4):
        if resolution > 3:
            self.resolution = resolution
        else:
            self.resolution = 3

        self.make_spherical()
        self.alpha      = alpha
        # generate plot
        self.figure          = plt.figure()
        self.subplot         = self.figure.add_subplot(111, projection='3d')
        self.subplot.set_anchor((0.75,0.75))

    def add_sphere(self, position = np.zeros(3),
                         radius  = 1.0,
                         **kwargs):
        if 'color' not in kwargs:
            kwargs['color'] = mplColors.rgb2hex(0.8*np.random.rand(3)) 
        if 'alpha' not in kwargs:
            kwargs['alpha'] = 0.2

        self.subplot.plot_surface(*(self.make_points(position,radius)), **kwargs)

    def add_polygon(self,points,**kwargs):
        if 'color' not in kwargs:
            kwargs['color'] = mplColors.rgb2hex(0.8*np.random.rand(3))

        polygon = Poly3DCollection([points], **kwargs)
        polygon.set_color(kwargs['color'])
        self.subplot.add_collection3d(polygon)

    def set_aspect(self,aspect_data):
        scale_from_ascpect = np.repeat([[np.min(aspect_data), np.max(aspect_data)]],3,axis=0)
        self.subplot.auto_scale_xyz(*scale_from_ascpect)

    @staticmethod
    def show(name=None):
        if name is None:
            plt.show()
        else:
            plt.savefig(name)
