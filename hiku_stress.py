#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 14:53:16 2019

Hiku stress calculations

@author: amt
"""

import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy import interpolate
from tde import tde
import utm

# read ascii mesh file
asc_file='/Users/amt/Documents/nz_repeaters/hiku.asc'
meshinfo=np.loadtxt(asc_file,skiprows=1)
centroidx, centroidy,_,_=utm.from_latlon(meshinfo[:,2], meshinfo[:,1])
node1x, node1y,_,_=utm.from_latlon(meshinfo[:,5], meshinfo[:,4])
node2x, node2y,_,_=utm.from_latlon(meshinfo[:,8], meshinfo[:,7])
node3x, node3y,_,_=utm.from_latlon(meshinfo[:,11], meshinfo[:,10])
centroidz=meshinfo[:,3]*1000 # all in meters
node1z=meshinfo[:,6]*1000
node2z=meshinfo[:,9]*1000
node3z=meshinfo[:,12]*1000

# plot the mesh
plt.figure()
plt.scatter(centroidx,centroidy,s=25,c=centroidz, cmap='jet')
plt.axis('equal')
plt.colorbar()

# read slip model
def read_laura(filename: str) -> dict:
   with open(filename, "r") as f:
       lines = f.read().split("\n")
   out_dict = None
   for line in lines[1:]:
       if len(line) == 0:
           continue
       out_dict = convert_line(line, out_dict)
   # Convert to numpy arrays
   for key, value in out_dict.items():
       out_dict.update({key: np.array(value)})
   return out_dict

def convert_line(line: str, out_dict: dict = None) -> dict:
   line = line[13:].split()
   out_dict = out_dict or dict(
       longitude=[], latitude=[], ztop=[], zbot=[], width=[], strike=[],
       dip=[], rake=[], U_ss=[], U_ds=[], U_te=[], U_tot=[], area=[])
   long, lat, ztop, zbot, width, strike, dip, rake, U_ss, U_ds, U_te, U_tot, area = (
       float(line[0]), float(line[1]), float(line[2]), float(line[3]), float(line[4]),
       float(line[5]), float(line[6]), float(line[7]), float(line[8]), float(line[9]),
       float(line[10]), float(line[11]), float(line[-1]))
   out_dict["longitude"].append(long)
   out_dict["latitude"].append(lat)
   out_dict["ztop"].append(ztop)
   out_dict["zbot"].append(zbot)
   out_dict["width"].append(width)
   out_dict["strike"].append(strike)
   out_dict["dip"].append(dip)
   out_dict["rake"].append(rake)
   out_dict["U_ss"].append(U_ss)
   out_dict["U_ds"].append(U_ds)
   out_dict["U_te"].append(U_te)
   out_dict["U_tot"].append(U_tot)
   out_dict["area"].append(area)
   return out_dict 

def akirich(strike, dip, rake):
    '''
    % strike (degrees clockwise from north), dip, and
    % rake (0=ll, 90=thrust, -90=normal, 180=RL)
    % I flipped these from pg 108 in A and R because I want x=east, y=north,
    % and z=up'''
    cosd = lambda x : np.cos( np.deg2rad(x) )
    sind = lambda x : np.sin( np.deg2rad(x) )
    n=[ sind(dip)*cosd(strike), -sind(dip)*sind(strike), cosd(dip) ]
    s=[ cosd(rake)*sind(strike)-cosd(dip)*sind(rake)*cosd(strike),
        cosd(rake)*cosd(strike)+cosd(dip)*sind(rake)*sind(strike),
        sind(rake)*sind(dip) ]  
    return n, s

slip_file='/Users/amt/Documents/nz_repeaters/g13a_info.out'
out_dict=read_laura(slip_file)
out_dict['x'], out_dict['y'],_,_=utm.from_latlon(out_dict['latitude'], out_dict['longitude'])
out_dict['ztop']*=1000
out_dict['zbot']*=1000

inds=np.where(out_dict['U_tot']>0)
sliponly={key: value[inds] for key, value in out_dict.items()}

sliponly['z']=(sliponly['ztop']+sliponly['zbot'])/2
plt.scatter(sliponly['x'],sliponly['y'],s=25,c=sliponly['z'], cmap='Blues')

# interpolate slip model to our mesh
slipinterp=interpolate.griddata(np.array([out_dict['x'],out_dict['y']]).T,out_dict['U_tot']/1000,
                        np.array([centroidx, centroidy]).T)
ssinterp=interpolate.griddata(np.array([out_dict['x'],out_dict['y']]).T,out_dict['U_ss']/1000,
                        np.array([centroidx, centroidy]).T)
dsinterp=interpolate.griddata(np.array([out_dict['x'],out_dict['y']]).T,out_dict['U_ds']/1000,
                        np.array([centroidx, centroidy]).T)
teinterp=interpolate.griddata(np.array([out_dict['x'],out_dict['y']]).T,out_dict['U_te']/1000,
                        np.array([centroidx, centroidy]).T)

# plot the mesh
plt.figure()
plt.scatter(centroidx, centroidy,s=25,c=dsinterp, cmap='jet')
plt.axis('equal')
plt.colorbar()

# define repeater location and properties
relat=-38.3064
relon=178.5861
redep=-13000
rex,rey,_,_=utm.from_latlon(relat, relon)
strike=33
dip=90
rake=-179
n,s=akirich(strike,dip,rake)
plt.plot(rex,rey,'ko')

# Calculate lame parameters
def lame_params(type1, val1, type2, val2):
    '''% Enter two elastic parameters of E=youngs modulus, G=shear modulus,
    % nu=poisson ratio, lambda, K=bulk modulus'''
    G=val1
    if type2=='nu':
        nu=val2;
        K=2*G*(1+nu)/(3*(1-2*nu))
        E=2*G*(1+nu)
        lamda=2*G*nu/(1-2*nu)
    elif type2=='lamda':
        lamda=val2
        K=lamda+2*G/3
        nu=lamda/(2*(lamda+G))
        E=G*(3*lamda+2*G)/(lamda+G)
    lame={}
    lame['K']=K
    lame['G']=G
    lame['E']=E
    lame['nu']=nu
    lame['lamda']=lamda
    return lame

lame=lame_params('G',30e9,'nu',0.25)
# Calculate stresses 
E=dict(xx=0.0, yy=0.0, zz=0.0, xy=0.0, xz=0.0, yz=0.0)

inds=np.where(slipinterp>0)
for ii in inds[0]:
    print(ii)
    tmp = tde.calc_tri_strains([rex], [rey], [redep], (node1x[ii],node1y[ii],node1z[ii]), (node2x[ii],node2y[ii],node2z[ii]), 
                               (node3x[ii],node3y[ii],node3z[ii]), lame['nu'], ssinterp[ii], teinterp[ii], dsinterp[ii])
    for key in E.keys():
        E[key] += tmp[key]
    print(tmp)
    print(E)
S = tde.strain_to_stress(E, lame['lamda'], lame['G'])