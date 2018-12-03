#!/usr/bin/python

from matplotlib import colors
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3D
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Voronoi, ConvexHull 
from mpl_toolkits.mplot3d import proj3d

def add_sphere(ax,positon = np.zeros(3), radius = 1.0, mycolor = colors.rgb2hex(np.random.rand(3)), **kwargs):
    u = np.linspace(0, 2 * np.pi, 20)
    v = np.linspace(0, np.pi, 10)
    x = positon[0] + 1* radius * np.outer(np.cos(u), np.sin(v))
    y = positon[1] + 1* radius * np.outer(np.sin(u), np.sin(v))
    z = positon[2] + 1* radius * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z, color=mycolor, alpha=0.4)

#basis = np.array([ [2, 0,0],
#                   [1,np.sqrt(3.),0],
#                   [0,0,4*np.sqrt(2.0/3.0)]])
#superCell = np.array([ [-0.5,-0.5,-0.5],
#                       [ 0.5, 0.5, 0.5]])
#superCell = np.array([ [ 0.0, 0.0, 0.0],
#                       [ 4.0/3.0, 2.0/3.0, 2.0*np.sqrt(2.0/3.0)]])

basis = np.array([
 [ 3.87753641405488 ,-0.00000000000000 ,-0.00000000000000],
 [-0.00000000000000 , 3.87753641405488 , 0.00000000000000],
 [ 0.00000000000000 , 0.00000000000000 ,12.09549987244948]
   ])

superCell = np.array([
 [ 0.50000000000000*3.87753641405488,  0.50000000000000*3.87753641405488 ,0.19398936278918*12.09549987244948],    
 [ 0.50000000000000*3.87753641405488,  0.50000000000000*3.87753641405488 ,0.80601063721082*12.09549987244948],    
 [ 0.00000000000000*3.87753641405488,  0.00000000000000*3.87753641405488 ,0.00000000000000*12.09549987244948],    
 [-0.00000000000000*3.87753641405488,  0.00000000000000*3.87753641405488 ,0.36270017366243*12.09549987244948],    
 [ 0.00000000000000*3.87753641405488, -0.00000000000000*3.87753641405488 ,0.63729982633757*12.09549987244948],    
 [ 0.00000000000000*3.87753641405488, -0.00000000000000*3.87753641405488 ,0.15023942974438*12.09549987244948],    
 [-0.00000000000000*3.87753641405488,  0.00000000000000*3.87753641405488 ,0.84976057025562*12.09549987244948],    
 [ 0.00000000000000*3.87753641405488,  0.50000000000000*3.87753641405488 ,0.37983766131592*12.09549987244948],    
 [ 0.50000000000000*3.87753641405488,  0.00000000000000*3.87753641405488 ,0.37983766131592*12.09549987244948],    
 [-0.00000000000000*3.87753641405488,  0.50000000000000*3.87753641405488 ,0.62016233868408*12.09549987244948],    
 [ 0.50000000000000*3.87753641405488, -0.00000000000000*3.87753641405488 ,0.62016233868408*12.09549987244948],    
 [ 0.50000000000000*3.87753641405488,  0.50000000000000*3.87753641405488 ,0.50000000000000*12.09549987244948]    
   ])
atomNames = [ 'Ba', 'Ba', 'Cu', 'Cu', 'Cu', 'O', 'O', 'O', 'O', 'O', 'O', 'Y' ]
center = np.zeros(3)
for a in superCell:
    center += a
center /= len(superCell)

cutOff = 13
multiplier = 1 + int(cutOff/np.min([np.max([np.abs(np.dot(e,d)) for e in basis]) for d in np.identity(3)]))

points = []
names  = []
for n,a in zip(atomNames,superCell):  
    points.append(a)
    names.append(n)
for i in range(-multiplier,multiplier+1):
    for j in range(-multiplier,multiplier+1):
        for k in range(-multiplier,multiplier+1):
            if i == 0 and j == 0 and k == 0:
                continue
            for ii,a in enumerate(superCell):  
                v = a + i*basis[0]+j*basis[1]+k*basis[2]
                if(np.linalg.norm(v-center) <= cutOff):
                    points.append(v)  
                    names.append(atomNames[ii])

points = np.array(points)
                   
# calculate spherical Voronoi diagram
sv = Voronoi(points)

# generate plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# plot generator points
for element in set(atomNames):
    random_color = colors.rgb2hex(np.random.rand(3))
    for name,atom in zip(names,points):
        if element == name:
            ax.scatter(atom[0], atom[1], atom[2], c=random_color, alpha=0.05)

# plot Voronoi vertices
#ax.scatter(sv.vertices[:, 0], sv.vertices[:, 1], sv.vertices[:, 2], c='g')


# indicate Voronoi regions (as Euclidean polygons)
regions = []
for i,atom in enumerate(superCell):
    regions.append(sv.regions[sv.point_region[i]])
for name,atom,region in zip(atomNames,superCell,regions):
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

    radius = cutOff*2;
    if len(vertices) > 3:
        convexHull = ConvexHull(vertices)
        random_color = colors.rgb2hex(np.random.rand(3))
        volume = convexHull.volume
        radiusA = np.power(3*volume/(4*np.pi),1.0/3.0)
        for simplex in convexHull.simplices:
            normal = np.cross(vertices[simplex[1]]-vertices[simplex[0]],vertices[simplex[2]]-vertices[simplex[0]])
            normal /= np.linalg.norm(normal)
            distance = np.abs(np.dot(normal,atom-vertices[simplex[0]]))
            if distance < radius:
                radius = distance
#            polygon = Poly3DCollection([vertices[simplex]], alpha=alpha)
#            polygon.set_color(random_color)
#            ax.add_collection3d(polygon)
        add_sphere(ax,atom,radius,random_color)   
        print name,radius,radiusA
    elif len(vertices) > 2:
        random_color = colors.rgb2hex(np.random.rand(3))
        polygon = Poly3DCollection([vertices], alpha=1.0)
        polygon.set_color(random_color)
        ax.add_collection3d(polygon)
    elif len(vertices) > 1:
        random_color = colors.rgb2hex(np.random.rand(3))
        line = Line3D(vertices[:,0],vertices[:,1],vertices[:,2], alpha=1.0, linewidth=0.1, color=random_color)
        ax.add_collection3d(line)
plt.show()
exit()
#            ax.plot_trisurf(sv.vertices[simplex][:, 0], sv.vertices[simplex][:, 1], sv.vertices[simplex][:, 2], color=colors.rgb2hex(np.random.rand(3)), linewidth=0, antialiased=False,  alpha=0.5)
angle = 0.0

while True:
    ax.view_init(np.sin(0.01*angle), angle)
    plt.draw()
    plt.pause(.001)
    angle += 3.5

#plt.show()
