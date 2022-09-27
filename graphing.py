# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 17:04:18 2022

@author: sulai
"""


import numpy as np
import pickle as pkl
import networkx as nx
import matplotlib.pyplot as plt



time=np.linspace(0,400,400)

params=[]
for i in range(5):
    
    
    p0=34+2*i
    
    for j in range(4):
        
        a0=9+j
        
        r_list=np.load("radius_array" + str(p0)+str(a0)+".npy")
        
        late_ratio=r_list[-1]/r_list[200]
        
        if late_ratio>1.01 or late_ratio<0.99:
            params.append([p0,a0])
        else:
            pass
        
        
        plt.plot(time,r_list)
        
        
print(params)