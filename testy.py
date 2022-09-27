# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 15:10:25 2022

@author: sulai
"""

import numpy as np
import pickle as pkl

import networkx as nx
from scipy.spatial import distance

import matplotlib.pyplot as plt

from tissue_new import tissue_generation


from cell_properties import  set_initial_cell_normals, area_perimeter_calculate, set_centroids
import force_functions_new as ff

import matplotlib.pyplot as plt


path='C:/Users/sulai/OneDrive - Imperial College London/PhD/Simulation/files/'

V,C=nx.read_gpickle(path+'stage_new.pickle'),nx.read_gpickle(path+'stagec_new.pickle')

xc=0
x_array=[]
y_array=[]
for i in C:
    x_array.append(xc)
    area=C.nodes[i]["area"]
    y_array.append(area)
    xc+=1
    print(area)
plt.scatter(x_array,y_array,c="r",marker="x")
#area_perimeter_calculate(V,C)

