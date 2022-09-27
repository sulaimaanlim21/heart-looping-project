# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 18:05:53 2021

@author: sulai
"""


import numpy as np
import pickle as pkl
import networkx as nx
from scipy.spatial import distance
import plotly.graph_objects as go
import plotly.express as px
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.pyplot import pause
import pylab



import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation








#vi=[62, 63, 64, 78, 79, 80, 166, 167, 168, 170, 171, 172, 182, 183, 184, 186, 187, 188]

step_number=0



for k in range(10):
    
    
    file_name = 't' + str(k)
    V = nx.read_gpickle(file_name + '.pickle')
    pos=nx.get_node_attributes(V,"pos")
    
    node_xyz = np.array([pos[v] for v in sorted(V)])
    edge_xyz = np.array([(pos[u], pos[v]) for u, v in V.edges()])
    
    # Create the 3D figure
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    
    # Plot the nodes - alpha is scaled by "depth" automatically
    ax.scatter(*node_xyz.T, s=30, ec="w")
    
    # Plot the edges
    for vizedge in edge_xyz:
        ax.plot(*vizedge.T, color="tab:gray")
    
    
    # def _format_axes(ax):
    #     """Visualization options for the 3D axes."""
    #     # Turn gridlines off
    #     ax.grid(False)
    #     # Suppress tick labels
    #     for dim in (ax.xaxis, ax.yaxis, ax.zaxis):
    #         dim.set_ticks([])
    #     # Set axes labels
    #     ax.set_xlabel("x")
    #     ax.set_ylabel("y")
    #     ax.set_zlabel("z")
    
    # ax.axis('auto')
    # _format_axes(ax)
    fig.tight_layout()
    fig.show()
         
    
     
    
    
    
   
    
    
    
    
    
    
    
    
    
    # trace_edges = go.Scatter3d(x=x_edges,
    #                         y=y_edges,
    #                         z=z_edges,
    #                         mode='lines',
    #                         line=dict(color='black', width=2))
    
    # trace_nodes = go.Scatter3d(x=x_nodes,
    #                          y=y_nodes,
    #                         z=z_nodes,
    #                         mode="markers",
    #                         marker=dict(size=2,color="blue"))
    
    
    
    # data = [trace_edges, trace_nodes]
    # fig = go.Figure(data=data)
    
    
    #fig.show()
    #fig.write_html('image' + str(k) +'.html', auto_open=True)
    
    
    
    
    #fig.write_image('t' + str(k)+'.png')
    plt.savefig(file_name + '.png')
    
    plt.clf()
    step_number+=1
    print(step_number)
    