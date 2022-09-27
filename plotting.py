# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 15:20:51 2022

@author: sulai
"""

import networkx as nx
import plotly.graph_objects as go
import plotly.express as px


def plot(V):
    pos=nx.get_node_attributes(V,"pos")
    
    length=len(pos)
    
    x_nodes = [pos[i][0] for i in V]# x-coordinates of nodes
    y_nodes = [pos[i][1] for i in V]# y-coordinates
    z_nodes = [pos[i][2] for i in V]# z-coordinates
    
    print("set positons")
    
    
    
    
    edge_list=V.edges()
    x_edges=[]
    y_edges=[]
    z_edges=[]
    
    #need to fill these with all of the coordiates
    for edge in edge_list:
        #format: [beginning,ending,None]
        x_coords = [pos[edge[0]][0],pos[edge[1]][0],None]
        
    
        y_coords = [pos[edge[0]][1],pos[edge[1]][1],None]
        
    
        z_coords = [pos[edge[0]][2],pos[edge[1]][2],None]
        
        x_edges += x_coords
        y_edges += y_coords
        z_edges += z_coords
     
        
    print("set edges")
    trace_edges = go.Scatter3d(x=x_edges,
                            y=y_edges,
                            z=z_edges,
                            mode='lines',
                            line=dict(color='black', width=2))
    
    trace_nodes = go.Scatter3d(x=x_nodes,
                             y=y_nodes,
                            z=z_nodes,
                            mode="markers",
                            marker=dict(size=2,color="blue"))
    
    
    
    data = [trace_edges, trace_nodes]
    fig = go.Figure(data=data)
    
    
    fig.show()
    fig.write_html('image.html', auto_open=True)