#!/usr/bin/python

from matplotlib import colors
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3D
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Voronoi, ConvexHull 
from mpl_toolkits.mplot3d import proj3d

# set input data

basis = np.array([ [2, 0,0],
                   [1,np.sqrt(3.),0],
                   [0,0,4*np.sqrt(2.0/3.0)]])

#superCell = np.array([ [-0.5,-0.5,-0.5],
#                       [ 0.5, 0.5, 0.5]])
superCell = np.array([ [ 0.0, 0.0, 0.0],
                       [ 4.0/3.0, 2.0/3.0, 2.0*np.sqrt(2.0/3.0)]])
center = np.zeros(3)
for a in superCell:
    center += a
center /= len(superCell)

cutOff = 5
multiplier = 1 + int(cutOff/np.min([np.max([np.abs(np.dot(e,d)) for e in basis]) for d in np.identity(3)]))

points = []
for a in superCell:  
    points.append(a)
for i in range(-multiplier,multiplier+1):
    for j in range(-multiplier,multiplier+1):
        for k in range(-multiplier,multiplier+1):
            if i == 0 and j == 0 and k == 0:
                continue
            for a in superCell:  
                v = a + i*basis[0]+j*basis[1]+k*basis[2]
                if(np.linalg.norm(v-center) <= cutOff):
                    points.append(v)  

points = np.array(points)
                   
# calculate spherical Voronoi diagram
sv = Voronoi(points)

# generate plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# plot generator points
ax.scatter(points[:, 0], points[:, 1], points[:, 2], c='b')

# plot Voronoi vertices
#ax.scatter(sv.vertices[:, 0], sv.vertices[:, 1], sv.vertices[:, 2], c='g')


# indicate Voronoi regions (as Euclidean polygons)
regions = []
for i,atom in enumerate(superCell):
    regions.append(sv.regions[sv.point_region[i]])
for region in regions:
    alpha = 0.4
    vertices = []
    average = np.zeros(3)
    ISINFINITY = False
    for point in region:
        if point != -1:
            vertices.append(sv.vertices[point])
            average += sv.vertices[point]
        else:
            ISINFINITY = True
#    if len(vertices) < 18:
#        alpha = 0.0
#    else:
#        alpha = 0.4
    if ISINFINITY:
        continue
        alpha = 0.05
        average /= len(vertices)
        vertices.append(average + 27720*(average-center))
    vertices = np.array(vertices)    
            
    if len(vertices) > 3:
        convexHull = ConvexHull(vertices)
        random_color = colors.rgb2hex(np.random.rand(3))
        for simplex in convexHull.simplices:
            polygon = Poly3DCollection([vertices[simplex]], alpha=alpha)
            polygon.set_color(random_color)
            ax.add_collection3d(polygon)
    elif len(vertices) > 2:
        random_color = colors.rgb2hex(np.random.rand(3))
        polygon = Poly3DCollection([vertices], alpha=1.0)
        polygon.set_color(random_color)
        ax.add_collection3d(polygon)
    elif len(vertices) > 1:
        random_color = colors.rgb2hex(np.random.rand(3))
        line = Line3D(vertices[:,0],vertices[:,1],vertices[:,2], alpha=1.0, linewidth=0.1, color=random_color)
        ax.add_collection3d(line)

#            ax.plot_trisurf(sv.vertices[simplex][:, 0], sv.vertices[simplex][:, 1], sv.vertices[simplex][:, 2], color=colors.rgb2hex(np.random.rand(3)), linewidth=0, antialiased=False,  alpha=0.5)
plt.show()
