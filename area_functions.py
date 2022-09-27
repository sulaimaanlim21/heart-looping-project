# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 12:14:03 2022

@author: sulai
"""

import numpy as np
import pickle as pkl

import networkx as nx
from scipy.spatial import distance

import matplotlib.pyplot as plt

from tissue_new import tissue_generation

from cell_properties import posi

import cell_properties as cp



def area_properties(V,C,c):
    vertex_list=cp.order_rh(V, C, c)
    
    area_total=np.array([0,0,0])
    n=len(vertex_list)
    
    for i in range(n):
        
        a_i=i
        b_i=(i+1) % n
        
        
    
        
        a=vertex_list[a_i]
        b=vertex_list[b_i]
        v_area=0.5*np.cross(cp.posi(V,a),cp.posi(V,b))
        
        area_total=area_total+v_area
        
    
    area_mag=distance.euclidean([0,0,0],area_total)
    area_vec=area_total
    
    unit_normal=(1/area_mag)*area_vec
    return area_mag,area_vec,unit_normal




def area_deriv_total(V,C,v):
    cell_list=V.nodes[v]["cells"]
    area_deriv_total=np.array([0,0,0])
    for c in cell_list:
        area_deriv=area_deriv_c(V, C, v, c)
        area_deriv_total=area_deriv_total+area_deriv
        
    return area_deriv_total
        

def area_deriv_c(V,C,v,c):
    ## area deriv is (1/A)*sum(Ak *dAk/dx)
    ## dAk/dx = (1/2)*(x_[i-1]-x_[i+1]) X ek    (ek is unit vector in k)
    ##call x+1 index b x-1 index a
    
    ex,ey,ez=np.array([1,0,0]),np.array([0,1,0]),np.array([0,0,1])
    A_mag,A_vec=area_test(V,C,c)
    vertex_list=cp.order_rh(V, C, c)
    n=len(vertex_list)
    i=vertex_list.index(v)
    b=vertex_list[(i+1)%n]
    a=vertex_list[i-1]
    
    
    
    
    
    dAx=(1/2)*A_vec[0]*(np.cross(np.array(posi(V,b))-np.array(posi(V,a)),ex))
    dAy=(1/2)*A_vec[1]*(np.cross(np.array(posi(V,b))-np.array(posi(V,a)),ey))
    dAz=(1/2)*A_vec[2]*(np.cross(np.array(posi(V,b))-np.array(posi(V,a)),ez))
    
    dA_by_dx=(1/A_mag)*(dAx+dAy+dAz)
    
    return np.array(dA_by_dx)



V,C=tissue_generation(16, 16)
cp.set_initial_cell_properties(V, C)