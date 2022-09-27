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

from force_functions2 import total_force_on_vertex, pressure_force2, bending_force2, slm_force, pressure_force, bending_force3, slm_force_anchor
from cell_properties import posi, get_centroid


dt=0.01

def sim_stage(V,C,parameters,name):
    ### parameters should look like [A1,A2,A3]



        
    for i in V:
        if V.nodes[i]["anchor"]==True:
            px,py,pz=posi(V,i)
    #         x_array=[]
    #         y_array=[]
    #         for cell in C:
    #             p1=get_centroid(V, C, cell)
    #             if pz-1.5*a<p1[2]<pz+1.5*a:
    #                 x_array.append(p1[0])
    #                 y_array.append(p1[1])
    #             else:
    #                 pass
                          
    #         xc,yc=np.mean(x_array),np.mean(y_array)
    #             #get angle
    #         x,y,z=posi(V,i)[0]-xc,posi(V,i)[1]-yc,posi(V,i)[2]
    #         radius=np.sqrt(x**2+y**2)
    #         ### set angle correctly
    #         if y>0:
    #             angle=np.arccos(x/radius)
    #         else:
    #             angle=2*np.pi-np.arccos(x/radius)
        
        
    #         x1,y1=0.7*radius*np.cos(angle)+xc,0.7*radius*np.sin(angle)+yc
    #         #V.nodes[i]["pos"]=np.array([x1,y1,z])
            
            V.nodes[i]["anchor_pos"]=np.array([px,py,pz])
            
    #     else:
    #         pass
    
        
    
    for k in range(100):
        pos = nx.get_node_attributes(V, "pos")
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
    
            else:
                pass
        #file_name = 't' + str(step_number)
    
        force_list = []
    
     
    
        for i in V:
          
                if V.nodes[i]["anchor"]==False:
                    force =total_force_on_vertex(V,C,i,parameters)+ bending_force3(V, C,  i, 3)+bending_force2(V, C, i, 3)
                else:
                    
                    
                    force = np.array([0,0,0])+total_force_on_vertex(V,C,i,parameters)+slm_force_anchor(V, i, 1)+bending_force3(V, C,  i, 3)
                    
                        
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
        
    
    
    nx.write_gpickle(V,'C:/Users/sulai/OneDrive - Imperial College London/PhD/Simulation/files/stage_' + str(name)+'.pickle')

    nx.write_gpickle(C,'C:/Users/sulai/OneDrive - Imperial College London/PhD/Simulation/files/stagec_' + str(name)+'.pickle')

    
    return (V,C)
