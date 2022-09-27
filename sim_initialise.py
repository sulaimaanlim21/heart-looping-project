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

from tissue_new import tissue_generation


from t1_transition import t1_transition
from asymmetric_growth import set_asymmetric_growth
from force_functions_new import total_force_on_vertex,bending_force2,grow_angle,bending_force3,grow_factor_z,bending_force_new
from cell_properties import posi, get_centroid
import cell_properties as cp
import force_functions_new as ff
from plotting import plot

KA=1
KP=1
eta=1
dt=0.01
#lt=0.015
A0=1
A1=0.9
P1=3.6*np.sqrt(A1)
a=2**(1/2)*3**(-3/4) 
P0=3.6
A2=1.5
P2=3.6*np.sqrt(A2)

step_number=0
sf=1
k1=1
k2=1


def sim_initialise(shape,tilt, name):

    n_x=shape[0]
    ny=shape[1]
    V,C=tissue_generation(n_x,ny)
    mid=ny*a*(3/4)
    lb=mid-2*a
    ub=mid+2*a
    max_y=ny*a*3/2
    ## add in the shrink
    
    cp.set_initial_cell_properties(V,C)
    
    for i in C:
        grow_angle(V,C,i,mid)
        grow_factor_z(V,C,i,mid)
        # for v in C.nodes[i]["vertices"]:
        #     set_asymmetric_growth(V, C, v, i)
            
     
    if tilt==True:
        for i in V:
            x,y,z=posi(V,i)
            if z>mid:
                y_shift=np.sin(((z-mid)/mid)*np.pi/2)*3
                V.nodes[i]["pos"]=[x,y+y_shift,z]
            else:
                y_shift=-np.sin(((mid-z)/mid)*np.pi/2)*3
                V.nodes[i]["pos"]=[x,y+y_shift,z]
    else:
        pass
    
    
    
    
    
    
    
    for i in V:
        V.nodes[i]["checked"]=False
    
    for i in C:
        cell_centre=get_centroid(V,C,i)
        vertices=C.nodes[i]["vertices"]
        if lb<cell_centre[2]<ub:
            
            C.nodes[i]["type"]="avc"
            for vert in vertices:
                V.nodes[vert]["type"]="avc"
                V.nodes[vert]["checked"]=True
        else:
            
            
            
            
            if cell_centre[2]<lb:
                C.nodes[i]["type"]="atr"
                for vert in vertices:
                    if V.nodes[vert]["checked"]==True:
                        pass
                    else:
                        V.nodes[vert]["type"]="atr"
            else:
                C.nodes[i]["type"]="ven"
                for vert in vertices:
                    if V.nodes[vert]["checked"]==True:
                        pass
                    else:
                        V.nodes[vert]["type"]="ven"
    
    
    
    
    for i in C:
        cell_centre=get_centroid(V,C,i)
        
        if lb<cell_centre[2]<ub:
            
            C.nodes[i]["type"]="avc"
            
        else:
        
            
            
            
            if cell_centre[2]<lb:
                C.nodes[i]["type"]="atr"
                
            else:
                C.nodes[i]["type"]="ven" 
    
        #grow_angle(V,C,i,mid)
              
    anchor_list=[]            
    for i in C:
        vertices=C.nodes[i]["vertices"]
        
        for v in vertices:
            if V.nodes[v]["anchor"]==True:
                anchor_list.append(i)
                
            else:
                pass
            
            
    for i in V:
        V.nodes[i]["drag"]=False
        
    for cell in anchor_list:
        verts=C.nodes[cell]["vertices"]
        for v in verts:
            V.nodes[v]["drag"]=True
    
  
    ##shrink top and bottom
    #for i in V:
        # if V.nodes[i]["anchor"]==True:
        #     px,py,pz=posi(V,i)
        #     x_array=[]
        #     y_array=[]
        #     for cell in C:
        #         p1=get_centroid(V, C, cell)
        #         if pz-1.5*a<p1[2]<pz+1.5*a:
        #             x_array.append(p1[0])
        #             y_array.append(p1[1])
        #         else:
        #             pass
                          
        #     xc,yc=np.mean(x_array),np.mean(y_array)
        #         #get angle
        #     x,y,z=posi(V,i)[0]-xc,posi(V,i)[1]-yc,posi(V,i)[2]
        #     radius=np.sqrt(x**2+y**2)
        #     ### set angle correctly
        #     if y>0:
        #         angle=np.arccos(x/radius)
        #     else:
        #         angle=2*np.pi-np.arccos(x/radius)
        
        
        #     x1,y1=0.7*radius*np.cos(angle)+xc,0.7*radius*np.sin(angle)+yc
        #     V.nodes[i]["pos"]=np.array([x1,y1,z])
            
            
            
        # else:
        #     pass
    
    
    for k in range(200):
        
     
             
        force_list=[]
        
     
        
        
        
    
        
        
        
        
        
            
        for i in V:
            
            
            
            
            
            # if V.nodes[i]["anchor"]==True:
            #     if posi(V,i)[2]>(3/4)*16*a:
            #         stretch_force=sf*np.array([0,0,1])
            #     else:
            #         stretch_force=-sf*np.array([0,0,1])
            # else:
            #     stretch_force=np.array([0,0,0])
                
            ## check if shrink cell
            if V.nodes[i]["anchor"]==False:
                force=bending_force3(V,C,i,1)
            else:
                force=bending_force3(V,C,i,1)
            #force=bending_force2(V,C,  i, 1)
            force_list.append([i,force])
            
            
            
            
            
            #print(force, "for node", i)
            
        
        for i in range(len(V)):
            vi=force_list[i][0]
            fi=np.array(force_list[i][1])
            
            current_pos=posi(V,vi)
            if V.nodes[vi]["anchor"]==False:
                new_pos=np.array(current_pos)+dt*np.array(fi)
            
            
            
            
            elif V.nodes[vi]["anchor"]==True:
                
                
                # if posi(V,vi)[2]<mid:
                
                      
                #     new_pos=current_pos#+dt*np.array([1/2,-1/2,0])
                    
                # else:
                     
                new_pos=current_pos#+dt*np.array([1/2,-1/2,0])
                
               
                
                #new_pos=np.array(current_pos)#+dt*np.array(fi)
            V.nodes[vi]["pos"]=new_pos
            
            
            
    
        
        
        
        
        
        
        
        
        
        
        
        # data = [trace_edges, trace_nodes]
        # fig = go.Figure(data=data)
        
        
        #fig.show()
        #fig.write_html('image' + str(step_number) +'.html', auto_open=True)
        
        
        #print(V)
        
        #fig.write_image('t' + str(k)+'.png')
        if k % 25 ==0:
            print(k)
        
        
        
    
    nx.write_gpickle(V,'C:/Users/sulai/OneDrive - Imperial College London/PhD/Simulation/files/init_' + str(name)+'.pickle')
    
    
    
    
    
    nx.write_gpickle(C,'C:/Users/sulai/OneDrive - Imperial College London/PhD/Simulation/files/initc_' +str(name)+'.pickle')
      
    return V,C
        
        
       