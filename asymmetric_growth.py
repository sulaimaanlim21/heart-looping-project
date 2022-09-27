# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 11:18:58 2022

@author: sulai
"""
import numpy as np
import pickle as pkl

import networkx as nx
from scipy.spatial import distance
import random
from shapely.geometry import Polygon
from tissue import tissue_generation
import matplotlib.pyplot as plt 

import cell_properties as cp
from cell_properties import posi



def set_asymmetric_growth(V,C,vertex_index,cell_index):
    vector=posi(V,vertex_index)-np.array(C.nodes[cell_index]["centroid"])
    mag=distance.euclidean(vector, [0,0,0])
    unit_vector=(1/mag)*vector
    
    dot=np.abs(np.dot(np.array([0,0,1]),unit_vector))
    
    V.nodes[vertex_index]["cell_"+str(cell_index)+"_growth"]=1+0.3*dot
    
    
    
    