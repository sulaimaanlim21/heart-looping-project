# -*- coding: utf-8 -*-
"""
Created on Tue Jan  4 16:22:41 2022

@author: sulai
"""

import numpy as np
import pickle as pkl

import networkx as nx
from scipy.spatial import distance

import matplotlib.pyplot as plt

from tissue import tissue_generation

from matplotlib.pyplot import pause
import pylab
from t1_transition import t1_transition
import matplotlib.animation as animation
from force_functions import total_force_on_vertex
from cell_properties import posi

P0=3.6
KA=1
KP=1
eta=1
dt=0.01
#lt=0.015
A0=1


V,C=tissue_generation(8,12)
step_number=0

for k in range(200):
    for i in V:
        neighbours=list(V[i])
        for j in neighbours:
            d_ij=distance.euclidean(posi(V,i),posi(V,j))
            
            if d_ij<0.1:
                
                t1_transition(V,C,i,j)
                
            else:
                pass
    #file_name = 't' + str(step_number) 
         
    force_list=[]
    
    # if k==5:
    #     for m in move_list:
        
    #         v_t=pos[m]
        
    #         v_move=[2,0]
    #         v_new=np.array(v_t)+np.array(v_move)
        
    #         V.nodes[m]["pos"]=v_new
    # else:
    #     pass
    ## updating position
    for i in V:
        

        force=total_force_on_vertex(V, C, i,P0,A0)
        force_list.append([i,force])
        
        
        #print(force, "for node", i)
        
    
    for i in range(len(V)):
        vi=force_list[i][0]
        fi=np.array(force_list[i][1])
        
        current_pos=posi(V,vi)
        if V.nodes[vi]["anchor"]==False:
            new_pos=np.array(current_pos)+dt*np.array(fi)
        else:
            new_pos=current_pos
        
        V.nodes[vi]["pos"]=new_pos
        
        
        
    for i in V:
        neighbours=list(V[i])
        for j in neighbours:
            d_ij=distance.euclidean(posi(V,i),posi(V,j))
            
            if d_ij<0.1:
                
                t1_transition(V,C,i,j)
                
            else:
                pass
        
    
   
    
    ## check for t1 transitions
    
    
    
    # print(t1_check)
    # print(t1_list)
    pos=nx.get_node_attributes(V,"pos")
    
    #plot=nx.draw_networkx(V,pos,node_size=20,node_color="blue",node_shape="+",with_labels=False)
    
    # for m in move_list:
    #     plt.scatter(pos[m][0],pos[m][1],c="r")
   
    
    print("time is", k)
    
    file_name = 't' + str(step_number) 
    nx.write_gpickle(V,file_name + '.pickle')
    
    
    
    
    file_name2 = 'tc' + str(step_number) 
    nx.write_gpickle(C,file_name2 + '.pickle')
    
    step_number+=1
