# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 18:17:07 2021

@author: sulai
"""



import numpy as np
import hexalattice
import pickle as pkl
from hexalattice.hexalattice import *


d=np.sqrt(2/np.sqrt(3))
a=2**(1/2)*3**(-3/4)
b=2**(-1/2)*3**(-1/4)
d0=[0,a]
d1=[b,a/2]
d2=[b,-a/2]
d3=[0,-a]
d4=[-b,-a/2]
d5=[-b,a/2]



def generate_cells(nx,ny):
    ## has to be multiple of 4
    cells=[]
    xs=d*nx/2
    ys=(ny/4)*2*a + (ny/4-1)*a + a/2
    
    hex_centers, _ = create_hex_grid(nx=nx,ny=ny,do_plot=False,min_diam=d)
    
    for i in range(nx*ny):
        cell={}
        
        shifted_center=[hex_centers[i][0]+xs,hex_centers[i][1]+ys]
        cell["center"]=shifted_center
        vertices=[]
        vertices.append([d0[0]+shifted_center[0],d0[1]+shifted_center[1]])
        vertices.append([d1[0]+shifted_center[0],d1[1]+shifted_center[1]])
        vertices.append([d2[0]+shifted_center[0],d2[1]+shifted_center[1]])
        vertices.append([d3[0]+shifted_center[0],d3[1]+shifted_center[1]])
        vertices.append([d4[0]+shifted_center[0],d4[1]+shifted_center[1]])
        vertices.append([d5[0]+shifted_center[0],d5[1]+shifted_center[1]])
        
        
        cell["vertices"]=vertices
    
    
        cells.append(cell)
    
    return cells

#cells=generate_cells(8)

#for i in range()
## amount to shift x by: d*n_x/2
## amount to shift y by: if n_y/2 is even:  ny=n/2 ny*2a + (ny-1)*a +a/2  