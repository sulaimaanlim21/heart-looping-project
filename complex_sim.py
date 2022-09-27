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

from cell_properties import posi, get_centroid
import plotly.graph_objects as go
import plotly.express as px



from force_functions import total_force_on_vertex,pressure_force,bending_force2
from cell_properties import posi, get_centroid


KA=1
KP=1
eta=1
dt=0.01
#lt=0.015
A0=1
A1=0.8
P1=3.4*np.sqrt(A1)
a=2**(1/2)*3**(-3/4) 
P0=3.4
A2=1.2
P2=3.4*np.sqrt(A2)

step_number=0
sf=1
k1=1
k2=1

## centre is (3/4)*ny*a
## shrink the target area of the central cells
## height is 2*a, overall height is ny*2*a, halfway point is ny*a, take say a cell either side, so range of 3a about centre? 

n_x=8
ny=16
V,C=tissue_generation(8,16)
mid=ny*a*(3/4)
lb=mid-2*a
ub=mid+5*a

## add in the shrink



for i in C:
    cell_centre=get_centroid(V,C,i)
    vertices=C.nodes[i]["vertices"]
    if lb<cell_centre[2]<ub:
        
        C.nodes[i]["type"]="avc"
        for vert in vertices:
            V.nodes[vert]["type"]="avc"
    else:
        
        
        
        
        if cell_centre[2]<lb:
            C.nodes[i]["type"]="atr"
            for vert in vertices:
                V.nodes[vert]["type"]="atr"
        else:
            C.nodes[i]["type"]="ven"
            for vert in vertices:
                V.nodes[vert]["type"]="ven"
            
# anchor_list=[]            
# for i in C:
#     vertices=C.nodes[i]["vertices"]
    
#     for v in vertices:
#         if V.nodes[v]["anchor"]==True:
#             anchor_list.append(i)
            
#         else:
#             pass
    
# for cell in anchor_list:
#     verts=C.nodes[cell]["vertices"]
#     for v in verts:
#         V.nodes[v]["anchor"]=True

### add in the twist in top half

# for i in V:
#     if V.nodes[i]["anchor"]==True and posi(V,i)[2]>12*a:
        
        
#         ##get angle
#         x,y,z=posi(V,i)[0],posi(V,i)[1],posi(V,i)[2]
#         radius=np.sqrt(x**2+y**2)
#         ### set angle correctly
#         if y>0:
#             angle=np.arccos(x/radius)
#         else:
#             angle=2*np.pi-np.arccos(x/radius)
            
#         new_angle=angle+np.pi/6
        
#         x1,y1=radius*np.cos(new_angle),radius*np.sin(new_angle)
#         V.nodes[i]["pos"]=np.array([x1,y1,z])
#     else:
#         pass


# ## add twist in bottom half
# for i in V:
#     if V.nodes[i]["anchor"]==True and posi(V,i)[2]<10*a:
        
        
#         ##get angle
#         x,y,z=posi(V,i)[0],posi(V,i)[1],posi(V,i)[2]
#         radius=np.sqrt(x**2+y**2)
#         ### set angle correctly
#         if y>0:
#             angle=np.arccos(x/radius)
#         else:
#             angle=2*np.pi-np.arccos(x/radius)
            
#         new_angle=angle-np.pi/6
        
#         x1,y1=radius*np.cos(new_angle),radius*np.sin(new_angle)
#         V.nodes[i]["pos"]=np.array([x1,y1,z])
#     else:
#         pass





for k in range(400):
    
    # A1=1-0.2*k/200
    # P1=6*a-0.4*k/200
    
    # A2=1+0.2*k/200
    # P2=6*a+0.4*k/200
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
    
    
    
    # for i in V:
    #     ## max height is 3/2 * ny * a
    #     #centre 1 is 15/2 * a + 3/4 * ny * a 
    #     #centre 2 is 9/2 * a 
        
    #     #direction is pos-centre
        
        
    #     c1=[0,0,25*a]
    #     c2=[0,0,5*a]
    #     if posi(V,i)[2]>20*a:
    #         df=np.array(posi(V,i))-np.array(c1)
    #         l=distance.euclidean(df,[0,0,0])
    #         d_f=np.array(df)/l
            
    #         pressure_force=d_f
    #     elif posi(V,i)[2]<10*a:
    #         df=np.array(posi(V,i))-np.array(c2)
    #         l=distance.euclidean(df,[0,0,0])
    #         d_f=np.array(df)/l
            
    #         pressure_force=d_f
    #     else:
    #         pressure_force=[0,0,0]
            
    
        
    ### pressure 
    ### get system centre
    # r=0.25
    # n_v=len(V)
    # pos=nx.get_node_attributes(V, "pos")
    # centre=np.array([0,0])
    # for i in V:
    #     centre=centre+(1/n_v)*np.array([posi(V,i)[0],posi(V,i)[1]])
        
    
    
    
    
    
    
        
    for i in V:
        
        
        
        
        
        # if V.nodes[i]["anchor"]==True:
        #     if posi(V,i)[2]>(3/4)*16*a:
        #         stretch_force=sf*np.array([0,0,1])
        #     else:
        #         stretch_force=-sf*np.array([0,0,1])
        # else:
        #     stretch_force=np.array([0,0,0])
            
        ## check if shrink cell
        if V.nodes[i]["type"]=="avc":
            force=total_force_on_vertex(V, C, i,P0,A0)+bending_force2(V,C,  i, 1)#+pressure_force(V, i, 1)
        elif V.nodes[i]["type"]=="atr":
            force=total_force_on_vertex(V, C, i,P2,A2)+bending_force2(V,C,  i, 1)#+pressure_force(V,i,1)
            
        else:
            force=total_force_on_vertex(V, C, i,P0,A0)+bending_force2(V,C,  i, 1)
        force_list.append([i,force])
        
        
        
        
        
        #print(force, "for node", i)
        
    
    for i in range(len(V)):
        vi=force_list[i][0]
        fi=np.array(force_list[i][1])
        
        current_pos=posi(V,vi)
        if V.nodes[vi]["anchor"]==False:
            new_pos=np.array(current_pos)+dt*np.array(fi)
        
        
        
        
        elif V.nodes[vi]["anchor"]==True:
            
            if k<150:    
                new_pos=current_pos+dt*np.array([0,0,fi[2]])
            else:
                new_pos=current_pos+dt*np.array([fi[0],fi[1],0])
        
        V.nodes[vi]["pos"]=new_pos
        


### pressure force
    
    ## check for t1 transitions
    
    
    
    # print(t1_check)
    # print(t1_list)
    pos=nx.get_node_attributes(V,"pos")
    
    length=len(pos)
    
    x_nodes = [pos[i][0] for i in V if V.nodes[i]["type"]!="avc"]# x-coordinates of nodes
    y_nodes = [pos[i][1] for i in V if V.nodes[i]["type"]!="avc"]# y-coordinates
    z_nodes = [pos[i][2] for i in V if V.nodes[i]["type"]!="avc"]# z-coordinates
    
    x_nodes2 = [pos[i][0] for i in V if V.nodes[i]["type"]=="avc"]# x-coordinates of nodes
    y_nodes2 = [pos[i][1] for i in V if V.nodes[i]["type"]=="avc"]# y-coordinates
    z_nodes2 = [pos[i][2] for i in V if V.nodes[i]["type"]=="avc"]# z-coordinates
    
    
      
    # x_nodes2 = [pos[i][0] for i in vi]# x-coordinates of nodes
    # y_nodes2 = [pos[i][1] for i in vi]# y-coordinates
    # z_nodes2 = [pos[i][2] for i in vi]
    
    
    edge_list=V.edges()
    x_edges=[]
    y_edges=[]
    z_edges=[]
    
    x_edges2=[]
    y_edges2=[]
    z_edges2=[]
    
    #need to fill these with all of the coordiates
    for edge in edge_list:
        #format: [beginning,ending,None]
        x_coords = [pos[edge[0]][0],pos[edge[1]][0],None]
        
    
        y_coords = [pos[edge[0]][1],pos[edge[1]][1],None]
        
    
        z_coords = [pos[edge[0]][2],pos[edge[1]][2],None]
        
        if V.nodes[edge[0]]["type"]=="avc" and V.nodes[edge[1]]["type"]=="avc":
        
            x_edges2 += x_coords
            y_edges2 += y_coords
            z_edges2 += z_coords
        else:
            x_edges += x_coords
            y_edges += y_coords
            z_edges += z_coords
            
     
        
    
    trace_edges = go.Scatter3d(x=x_edges,
                            y=y_edges,
                            z=z_edges,
                            mode='lines',
                            line=dict(color='black', width=1))
    
    trace_edges2 = go.Scatter3d(x=x_edges2,
                            y=y_edges2,
                            z=z_edges2,
                            mode='lines',
                            line=dict(color='red', width=1))
    
    trace_nodes = go.Scatter3d(x=x_nodes,
                             y=y_nodes,
                            z=z_nodes,
                            mode="markers",
                            marker=dict(size=1.5,color="blue")
                            )
    
    trace_nodes2 = go.Scatter3d(x=x_nodes2,
                             y=y_nodes2,
                            z=z_nodes2,
                            mode="markers",
                            marker=dict(size=2,color="red")
                            )
    
    
    
    data = [trace_edges, trace_nodes,trace_nodes2,trace_edges2]
    fig = go.Figure(data=data)
    
    
    
    
    # layout = go.Layout(
    #     width=1024,
    #     height=1024,
    #     scene=dict(
    #                aspectmode='data', #this string can be 'data', 'cube', 'auto', 'manual'
    #                #a custom aspectratio is defined as follows:
    #                aspectratio=dict(x=1, y=1, z=1)
    #                )
    #     )
    
    
    
    fig['layout']['scene']['aspectmode']='manual'
    fig['layout']['scene']['aspectratio']=dict(x=1,y=1,z=1)
    
    # data = [trace_edges, trace_nodes]
    # fig = go.Figure(data=data)
    
    
    #fig.show()
    #fig.write_html('image' + str(k) +'.html', auto_open=True)
    
    
    
    
    # pio.write_image(fig, 't' + str(step_number)+'.png', format=None,
    #                   scale=None, width=None, height=None)
    
   
    fig.write_image('t' + str(step_number)+'.png',engine='orca') 
    print(k)
   
