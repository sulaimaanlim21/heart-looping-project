# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 11:10:08 2022

@author: sulai
"""

import numpy as np
import pickle as pkl

import networkx as nx
from scipy.spatial import distance

import matplotlib.pyplot as plt

from tissue import tissue_generation
import plotly.graph_objects as go
import plotly.express as px


from t1_transition import t1_transition

from force_functions_new import total_force_on_vertex, bending_force2, bending_force3, slm_force_anchor,bending_force_new
from cell_properties import posi, get_centroid
import cell_properties as cp
import force_functions_new as ff

dt=0.01
a=2**(1/2)*3**(-3/4)
def sim_bend(V,C,name):
    ### parameters should look like [A1,A2,A3]


    for i in V:
        if V.nodes[i]["anchor"] == False:
            neighbours = list(V[i])
            for j in neighbours:
                d_ij = distance.euclidean(posi(V, i), posi(V, j))
                if V.nodes[i]["type"] == "avc":
                    if d_ij < 0.01:
    
                        t1_transition(V, C, i, j, 0.01)
    
                    else:
                        pass
                else:
                    if d_ij < 0.05:
    
                        t1_transition(V, C, i, j, 0.05)
    
                    else:
                       pass         
        
    
    
        

    for k in range(50):
        bend_k=0
        avm_k=0
        cp.set_cell_properties(V, C)
        
        
        
        
        
        pos = nx.get_node_attributes(V, "pos")
       
        
    
        force_list = []
    
     
    
        for i in V:
          
                if V.nodes[i]["anchor"]==False:
                    
                    #avm_force=total_force_on_vertex(V, C, i, parameters)
                    bend_force=bending_force_new(V, C, i)#+bending_force3(V, C, i, 1)
                    force =bend_force
                    
                    
                    
                else:
                    
                    
                    force = np.array([0,0,0])
                    
                        
                force_list.append([i, force])
    
        for i in range(len(V)):
            vi = force_list[i][0]
            fi = np.array(force_list[i][1])
    
            current_pos = posi(V, vi)
            if V.nodes[vi]["anchor"] == False:
                new_pos = np.array(current_pos)+dt*np.array(fi)
    
            elif V.nodes[vi]["anchor"] == True:
               
    
               new_pos = current_pos+dt*np.array(fi)
                
                    
    
            V.nodes[vi]["pos"] = new_pos
    
            pos = nx.get_node_attributes(V, "pos")
            
        if k % 25 ==0:
            print(k)
        else:
            pass
    
        
        
        
    nx.write_gpickle(V,'C:/Users/sulai/OneDrive - Imperial College London/PhD/Simulation/files/stage_' + str(name)+'.pickle')

    nx.write_gpickle(C,'C:/Users/sulai/OneDrive - Imperial College London/PhD/Simulation/files/stagec_' + str(name)+'.pickle')

  
    return (V,C)
