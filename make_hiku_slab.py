from numpy import genfromtxt,where,float64,r_,expand_dims,c_,zeros
from glob import glob
import sys

sys.path.append('/Users/amt/Documents/nz_repeaters')
from gmsh_tools import xyz2gmsh
import matplotlib.pyplot as plt
import numpy as np

#Gmsh output file name
gmsh_out=u'/Users/amt/Documents/nz_repeaters/hiku.gmsh'

# Load contour file
contours=genfromtxt('/Users/amt/Documents/nz_repeaters/hiku_cont_wtrench.txt')

#Depth filter
maxdepth=-150
contours=contours[np.where(contours[:,2]>maxdepth)[0],:]

# Plot them
plt.figure()
plt.scatter(contours[:,0],contours[:,1],s=25,c=contours[:,2], cmap='jet')
plt.axis('equal')
plt.colorbar()

#Line filters
[L1x1 , L1y1] = [177.2,  -41]
[L1x2 , L1y2] = [174, -39]

#By regions
#Equations of the line filters
m1=(L1y1-L1y2)/(L1x1-L1x2)
b1=L1y1-m1*L1x1
#Get above line 1
ytest=m1*contours[:,0]+b1
i=where((ytest-contours[:,1])<=0)[0] 
contours=contours[i,:]


plt.figure()
plt.scatter(contours[:,0],contours[:,1],s=25,c=contours[:,2], cmap='jet')
plt.axis('equal')
plt.colorbar()

##Write gmsh file
#xyz2gmsh(gmsh_out,contours[:,0],contours[:,1],contours[:,2],coord_type='UTM',projection_zone='60H')