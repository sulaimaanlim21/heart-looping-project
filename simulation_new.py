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

from force_functions import total_force_on_vertex, pressure_force2, bending_force2, slm_force, pressure_force, bending_force3, slm_force_anchor
from cell_properties import posi, get_centroid


KA = 1
KP = 1
eta = 1
dt = 0.01
# lt=0.015
A0 = 1
A1 = 0.7
P1 = 3.6*np.sqrt(A1)
a = 2**(1/2)*3**(-3/4)
P0 = 3.6
A2 = 1.4
P2 = 3.6*np.sqrt(A2)
A3 = 2
P3 = 3.6*np.sqrt(A3)

step_number = 0
sf = 1
n_x = 16
ny = 16

mid = ny*a*(3/4)
lb = mid-3*a
ub = mid+3*a

n=199
V = nx.read_gpickle('t_esa4' + str(n) + '.pickle')
C = nx.read_gpickle('t_esa4c' + str(n) +'.pickle')


# for i in C:
#     grow_angle(V,C,i,mid)
        
for i in V:
    if V.nodes[i]["anchor"]==True:
        px,py,pz=posi(V,i)
        x_array=[]
        y_array=[]
        for cell in C:
            p1=get_centroid(V, C, cell)
            if pz-1.5*a<p1[2]<pz+1.5*a:
                x_array.append(p1[0])
                y_array.append(p1[1])
            else:
                pass
                      
        xc,yc=np.mean(x_array),np.mean(y_array)
            #get angle
        x,y,z=posi(V,i)[0]-xc,posi(V,i)[1]-yc,posi(V,i)[2]
        radius=np.sqrt(x**2+y**2)
        ### set angle correctly
        if y>0:
            angle=np.arccos(x/radius)
        else:
            angle=2*np.pi-np.arccos(x/radius)
    
    
        x1,y1=0.7*radius*np.cos(angle)+xc,0.7*radius*np.sin(angle)+yc
        #V.nodes[i]["pos"]=np.array([x1,y1,z])
        
        V.nodes[i]["anchor_pos"]=np.array([x,y,z])
        
    else:
        pass

    

for k in range(200):
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

    # for i in V:
    # if k>200:
    #     if V.nodes[i]["type"]=="avc" or V.nodes[i]["anchor"]==True:
    #         twist_force=np.array([0,0,0])
    #     else:
    #         pz=pos[i][2]
    #         yd=distance.euclidean(mid,pz)
    #         twist_factor=np.sin(np.pi*((yd+2*a)/mid))
    #         if twist_factor>0.2:
    #             x_array=[]
    #             y_array=[]
    #             for p1 in V:
    #                 if pz-a/2<pos[p1][2]<pz+a/2:
    #                     x_array.append(pos[p1][0])
    #                     y_array.append(pos[p1][0])
    #                 else:
    #                     pass

    #             xc,yc=np.mean(x_array),np.mean(y_array)

    #             if  posi(V,i)[2]>mid:

    #                 ##get angle
    #                 x,y,z=posi(V,i)[0]-xc,posi(V,i)[1]-yc,posi(V,i)[2]
    #                 radius=np.sqrt(x**2+y**2)
    #                 ### set angle correctly
    #                 if y>0:
    #                     angle=np.arccos(x/radius)
    #                 else:
    #                     angle=2*np.pi-np.arccos(x/radius)

    #                 new_angle=angle+twist_factor*np.pi/(600*3)

    #                 x1,y1=radius*np.cos(new_angle)+xc,radius*np.sin(new_angle)+yc

    #             else:
    #                 x,y,z=posi(V,i)[0]-xc,posi(V,i)[1]-yc,posi(V,i)[2]
    #                 radius=np.sqrt(x**2+y**2)
    #                 ### set angle correctly
    #                 if y>0:
    #                     angle=np.arccos(x/radius)
    #                 else:
    #                     angle=2*np.pi-np.arccos(x/radius)

    #                 new_angle=angle-twist_factor*np.pi/(600*3)

    #                 x1,y1=radius*np.cos(new_angle)+xc,radius*np.sin(new_angle)+yc

    #             V.nodes[i]["pos"]=[x1,y1,z]
    #         else:
    #             pass
    # else:
    #     pass

    for i in V:
      
            if V.nodes[i]["anchor"]==False:
                force =total_force_on_vertex(V,C,i)+ bending_force3(V, C,  i, 1)+bending_force2(V, C, i, 1)
            else:
                if V.nodes[i]["pos"][2]>mid:
                
                    force = np.array([0,0,0])+total_force_on_vertex(V,C,i)+slm_force_anchor(V, i, 1)+bending_force3(V, C,  i, 1)
                else:
                    force = np.array([0,0,0])+total_force_on_vertex(V,C,i)+slm_force_anchor(V, i, 1)+bending_force3(V, C,  i, 1)

            force_list.append([i, force])

    for i in range(len(V)):
        vi = force_list[i][0]
        fi = np.array(force_list[i][1])

        current_pos = posi(V, vi)
        if V.nodes[vi]["anchor"] == False:
            new_pos = np.array(current_pos)+dt*np.array(fi)

        elif V.nodes[vi]["anchor"] == True:
            if posi(V,vi)[2]>mid:

                new_pos = current_pos+dt*np.array(fi)
            else:
                new_pos = current_pos+dt*np.array(fi)

        V.nodes[vi]["pos"] = new_pos

        pos = nx.get_node_attributes(V, "pos")

    if k == 50 or k == 100 or k == 150 or k == 200 or k == 250 or k == 300 or k == 350:
        length = len(pos)

        for i in V:
            V.nodes[i]["checked"] = False

        for i in C:

            vertices = C.nodes[i]["vertices"]
            if C.nodes[i]["type"] == "avc":

                for vert in vertices:
                    V.nodes[vert]["type"] = "avc"
                    V.nodes[vert]["checked"] = True
            else:
                for vert in vertices:
                    if V.nodes[vert]["checked"] == True:
                        pass
                    else:
                        V.nodes[vert]["type"] = "not"
                        V.nodes[vert]["checked"] = True

        x_nodes = [pos[i][0] for i in V if V.nodes[i]
                   ["type"] != "avc"]  # x-coordinates of nodes
        y_nodes = [pos[i][1]
                   for i in V if V.nodes[i]["type"] != "avc"]  # y-coordinates
        z_nodes = [pos[i][2]
                   for i in V if V.nodes[i]["type"] != "avc"]  # z-coordinates

        x_nodes2 = [pos[i][0] for i in V if V.nodes[i]
                    ["type"] == "avc"]  # x-coordinates of nodes
        y_nodes2 = [pos[i][1]
                    for i in V if V.nodes[i]["type"] == "avc"]  # y-coordinates
        z_nodes2 = [pos[i][2]
                    for i in V if V.nodes[i]["type"] == "avc"]  # z-coordinates

        # x_nodes2 = [pos[i][0] for i in vi]# x-coordinates of nodes
        # y_nodes2 = [pos[i][1] for i in vi]# y-coordinates
        # z_nodes2 = [pos[i][2] for i in vi]

        edge_list = V.edges()
        x_edges = []
        y_edges = []
        z_edges = []

        x_edges2 = []
        y_edges2 = []
        z_edges2 = []

        # need to fill these with all of the coordiates
        for edge in edge_list:
            #format: [beginning,ending,None]
            x_coords = [pos[edge[0]][0], pos[edge[1]][0], None]

            y_coords = [pos[edge[0]][1], pos[edge[1]][1], None]

            z_coords = [pos[edge[0]][2], pos[edge[1]][2], None]

            if V.nodes[edge[0]]["type"] == "avc" and V.nodes[edge[1]]["type"] == "avc":

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
                                   line=dict(color='blue', width=1))

        trace_edges2 = go.Scatter3d(x=x_edges2,
                                    y=y_edges2,
                                    z=z_edges2,
                                    mode='lines',
                                    line=dict(color='black', width=1))

        trace_nodes = go.Scatter3d(x=x_nodes,
                                   y=y_nodes,
                                   z=z_nodes,
                                   mode="markers",
                                   marker=dict(size=1.5, color="blue")
                                   )

        trace_nodes2 = go.Scatter3d(x=x_nodes2,
                                    y=y_nodes2,
                                    z=z_nodes2,
                                    mode="markers",
                                    marker=dict(size=2, color="red")
                                    )

        data = [trace_edges, trace_edges2]
        fig = go.Figure(data=data)

        fig['layout']['scene']['aspectmode'] = 'data'

        # fig.show()
        fig.write_html('image' + str(k) + '.html', auto_open=True)

    file_name = 't_esa5' + str(n+k+1)
    nx.write_gpickle(V, file_name + '.pickle')

    file_name2 = 't_esa5c' + str(n+k+1)
    nx.write_gpickle(C, file_name2 + '.pickle')

    step_number += 1
    print(str(n+k+1))
