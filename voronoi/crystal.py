#!/usr/bin/python3

from matplotlib import colors
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3D
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Voronoi, ConvexHull 
from mpl_toolkits.mplot3d import proj3d
from sys import argv,path
path.insert(0,'../')
from JorG.loadsave import load_POSCAR

if len(argv) == 1:
    POSCARname = "POSCAR"
else:
    POSCARname = argv[1]

def add_sphere(ax,positon = np.zeros(3), radius = 1.0, color=colors.rgb2hex(0.8*np.random.rand(3)), res=3, alpha = 0.4, **kwargs):
    u = np.linspace(0, 2*np.pi, res)
    v = np.linspace(0, np.pi, res)
    x = positon[0] + radius * np.outer(np.cos(u), np.sin(v))
    y = positon[1] + radius * np.outer(np.sin(u), np.sin(v))
    z = positon[2] + radius * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z, color=color, alpha=alpha, **kwargs)

def support(center,radius,direction):
    project = np.dot(center,direction)
    return [project + radius, project - radius]

def supportXYZ(center, radius):
    differences = []
    for v in np.identity(3):
        differences.append(support(center,radius,v))
    differences = np.array(differences)    
    return [np.min(differences),np.max(differences)]

data = load_POSCAR(POSCARname)
basis = data['directions']
superCell = np.array([ x[1] for x in data['cell']])
atomNames = [ data['atomNames'][atom[0]] for atom in data['cell'] ]

atomColors = {}
for name in set(atomNames):
    atomColors[name] = colors.rgb2hex(0.8*np.random.rand(3))

center = np.mean(superCell, axis=0)

cutOff = np.linalg.norm(np.sum(basis,axis=0))

points = []
names  = []
for n,a in zip(atomNames,superCell):  
    points.append(a)
    names.append(n)

multipliers = [-1,0,1]
for i in multipliers:
    for j in multipliers:
        for k in multipliers:
            if i == 0 and j == 0 and k == 0:
                continue
            for name,a in zip(atomNames,superCell):  
                v = a + i*basis[0]+j*basis[1]+k*basis[2]
                if(np.linalg.norm(v-center) <= cutOff):
                    points.append(v)  
                    names.append(name)
points = np.array(points)
                   
# calculate Voronoi diagram
sv = Voronoi(points)

# generate plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_anchor((0.75,0.75))

# indicate Voronoi regions (as Euclidean polygons)
regions = []
for i,atom in enumerate(superCell):
    regions.append(sv.regions[sv.point_region[i]])

aspect_data = []
for name,atom,region in zip(atomNames,superCell,regions):
    alpha = 0.4
    vertices = []
    for point in region:
        if point != -1:
            vertices.append(sv.vertices[point])
    vertices = np.array(vertices)    

    radius = cutOff*2;
    if len(vertices) > 3:
        convexHull = ConvexHull(vertices)
        volume = convexHull.volume
        WSradius = np.power(3*volume/(4*np.pi),1.0/3.0)
        for simplex in convexHull.simplices:
            normal = np.cross(vertices[simplex[1]]-vertices[simplex[0]],vertices[simplex[2]]-vertices[simplex[0]])
            normal /= np.linalg.norm(normal)
            distance = np.abs(np.dot(normal,atom-vertices[simplex[0]]))
            if distance < radius:
                radius = distance
            polygon = Poly3DCollection([vertices[simplex]], alpha=0.2)
            polygon.set_color(atomColors[name])
            ax.add_collection3d(polygon)
        add_sphere(ax,atom,radius,atomColors[name],alpha=1.0,res=10)   
        aspect_data.append(supportXYZ(atom,radius))
        print(name,radius,WSradius)

scale_from_ascpect = [np.min(aspect_data), np.max(aspect_data)]
ax.auto_scale_xyz(scale_from_ascpect,scale_from_ascpect,scale_from_ascpect)        
plt.show()
exit()

