# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 18:05:53 2021

@author: sulai
"""


import numpy as np
import pickle as pkl
import networkx as nx
from scipy.spatial import distance
import plotly
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
#plotly.io.orca.config.executable = '/path/to/orca'
#plotly.io.orca.config.save()
#pio.kaleido.scope.mathjax = None
step_number=0

path='C:/Users/sulai/OneDrive - Imperial College London/PhD/Simulation/big_files/'
file_name = 'stage_new2'
file_name2 = 'stagec_new2'
    
for k in range(200):
    
    
    V = nx.read_gpickle(path+file_name + str(200+k)+'.pickle')
    pos=nx.get_node_attributes(V,"pos")
    
    length=len(pos)
    
    x_nodes = [pos[i][0] for i in V ]# x-coordinates of nodes
    y_nodes = [pos[i][1] for i in V ]# y-coordinates
    z_nodes = [pos[i][2] for i in V ]# z-coordinates
    
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
        
        # if V.nodes[edge[0]]["type"]=="avc" and V.nodes[edge[1]]["type"]=="avc":
        
            # x_edges2 += x_coords
            # y_edges2 += y_coords
            # z_edges2 += z_coords
        #else:
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
                            marker=dict(size=0.7,color="blue")
                            )
    
    trace_nodes2 = go.Scatter3d(x=x_nodes2,
                             y=y_nodes2,
                            z=z_nodes2,
                            mode="markers",
                            marker=dict(size=0.7,color="red")
                            )
    
    
    
    data = [trace_edges]
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
    
    
    
    
    camera = dict(
    eye=dict(x=2, y=0, z=0))

    fig['layout']['scene']['aspectmode']='data'
    fig.update_layout(scene_camera=camera)
    fig.update_layout(
    scene = dict(
        xaxis = dict(range=[-8,8],),
                     yaxis = dict( range=[-8,8],),
                     zaxis = dict( range=[-2,20],)))
    fig.update_xaxes(showticklabels=False) # Hide x axis ticks 
    fig.update_yaxes(showticklabels=False)
    # fig = go.Figure(data=data)
    
    
    #fig.show()
    #fig.write_html('image' + str(k) +'.html', auto_open=True)
    
    
    print(k)
    
    
    #pio.write_image(fig, 't' + str(step_number)+'.png', format=None,
                      #scale=None, width=None, height=None)
    
    
    fig.write_image(path+file_name + str(k)+'.png',engine='orca')
    #plt.savefig(file_name + '.png')
    
    #plt.clf()
      
    #print("image",k)
    step_number+=1