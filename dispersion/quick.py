#!/usr/bin/python

import numpy as np
from numpy.linalg import norm
import scipy as sp

RES    = 1000
a      = np.sqrt(3)
R      = np.sqrt(np.sqrt(2)+np.sqrt(5))
dMax   = np.power(a*1.73205,2)+R*R
eps    = 1e-5
bandNo = 2

# band: [bandid,[position]]
elementaryCell = [
    [0,np.array([0,0,0])],
    [1,np.array([0,0,R])]
    ]

latticeDir = np.array([
    [a,0,0],
#    [0,a,0]
    [0.5*a,np.sqrt(0.75)*a,0]
    ])

crystal   = []
crystalID = []

hoppingsList = [{
R : -1.054548660071587429e+00,
np.sqrt(np.power(a*1,2)) : -9.111423709522717407e-02,
np.sqrt(np.power(a*1,2)+R*R) : 1.202010422660697471e-02,
np.sqrt(np.power(a*2.64575,2)) : -1.856634136686042061e-03,
np.sqrt(np.power(a*2.64575,2)+R*R) : -6.536540256278850543e-04,
np.sqrt(np.power(a*2,2)) : 8.118475400165508366e-03,
np.sqrt(np.power(a*2,2)+R*R) : 8.670651434340008760e-04,
np.sqrt(np.power(a*1.73205,2)) : 1.127519964026840044e-02,
np.sqrt(np.power(a*1.73205,2)+R*R) : 7.243177158263955926e-04
},{
R : -9.263369982013420767e-01,
np.sqrt(np.power(a*1,2)) : -3.322024249096183768e-02,
np.sqrt(np.power(a*1,2)+R*R) : -9.769361240990257975e-03,
np.sqrt(np.power(a*2.64575,2)) : -5.795721794059642031e-03,
np.sqrt(np.power(a*2.64575,2)+R*R) : -1.387567976617446125e-03,
np.sqrt(np.power(a*2,2)) : 1.399067648380310398e-02,
np.sqrt(np.power(a*2,2)+R*R) : 2.400531866005962273e-03,
np.sqrt(np.power(a*1.73205,2)) : 1.669508295861291819e-02,
np.sqrt(np.power(a*1.73205,2)+R*R) : -5.636634357015664169e-03
},{
R : -1.266491784841059065e+00,
np.sqrt(np.power(a*1,2)	-1.099925452) : 22695618e-01,
np.sqrt(np.power(a*1,2)+R*R) : 6.828880293392308991e-04,
np.sqrt(np.power(a*2.64575,2)) : -1.076173667481636689e-02,
np.sqrt(np.power(a*2.64575,2)+R*R) : -2.533049168704091386e-03,
np.sqrt(np.power(a*2,2)) : 2.677449119804332228e-02,
np.sqrt(np.power(a*2,2)+R*R) : 7.665076185820022831e-04,
np.sqrt(np.power(a*1.73205,2)) : 3.910870464547556991e-02,
np.sqrt(np.power(a*1.73205,2)+R*R) : -7.598662885085487491e-03
},{
R : -3.748791100244428659e-02,
np.sqrt(np.power(a*1,2)) : -3.299327041564717455e-01,
np.sqrt(np.power(a*1,2)+R*R) : -6.361317995109867784e-03,
np.sqrt(np.power(a*2.64575,2)) : -2.234717286063514149e-02,
np.sqrt(np.power(a*2.64575,2)+R*R) : -5.502882493887413344e-04,
np.sqrt(np.power(a*2,2)) : 5.137770464547557203e-02,
np.sqrt(np.power(a*2,2)+R*R) : 1.126983667481636607e-03,
np.sqrt(np.power(a*1.73205,2)) : 8.163561320293188928e-02,
np.sqrt(np.power(a*1.73205,2)+R*R) : 1.750041491442503152e-03
}]


for i in range(-10,10+1):
  for j in range(-40,40+1):
    for atom in elementaryCell:
      dMin = np.inf # minimal distance to the elementary cell
      newAtom = [atom[0],atom[1] + i*latticeDir[0] + j*latticeDir[1]] # atom translation
      for atom2 in elementaryCell: # checking minimal distance
        d = norm(newAtom[1] - atom[1])
        if d < dMin:
          dMin = d  
      if dMin <= dMax + eps:
        crystal.append(newAtom)  

momenta = []
for kx in np.linspace(0,np.pi,RES):
    momenta.append(np.array([kx,0.0,0.0]))
for ky in np.linspace(0,np.pi,RES):
    momenta.append(np.array([np.pi,ky,0.0]))
for ks in np.linspace(np.pi,0,RES):
    momenta.append(np.array([ks,ks,0.0]))

for i,hoppings in zip(range(len(hoppingsList)),hoppingsList):
  with open("band%d.dat"%i,"w+") as f:
    maxE = np.zeros(bandNo)
    minE = np.zeros(bandNo)
    for j in range(bandNo):
        minE[j] = np.inf
        maxE[j] =-np.inf
    for momentum in momenta:
      H = 0.j*np.zeros((bandNo,bandNo))
      for band in elementaryCell:
          for atom in crystal:
              distance = norm(band[1] - atom[1])
              H[band[0],atom[0]] += np.exp(-1.j*np.dot(momentum,band[1] - atom[1])/a)*hoppings.get(distance, 0.0)
      sol = np.linalg.eigh(H) 
      for k in momentum:
          f.write("%.15f "%k)
      for j,E in enumerate(sol[0]):
          if E > maxE[j]:
              maxE[j] = E
          if E < minE[j]:
              minE[j] = E
          f.write("%.15f "%E)
      f.write("\n")
      
    for i,(eP,eM) in enumerate(zip(maxE,minE)):
        print i, eP, eM, eP - eM,
    print ""    
