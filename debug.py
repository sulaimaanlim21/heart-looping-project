# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 14:57:42 2022

@author: sulai
"""

import numpy as np
import pickle as pkl
import networkx as nx
from scipy.spatial import distance
import plotly.graph_objects as go
import plotly.express as px

from matplotlib.pyplot import pause
import pylab
from plotting import plot


import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation

import cell_properties as cp
import force_functions as ff
import matplotlib.pyplot as plt
from shapely.geometry import Polygon


P0=3.4



##### where i have got to: am checking why the z force is so big for these vertices, and the perimeter is bigger than it should be for the cells in cell list



V = nx.read_gpickle('t4.pickle')
C=nx.read_gpickle('tc4.pickle')
pos=nx.get_node_attributes(V,"pos")


vi,fi=[],[]
for i in V:
    force=ff.total_force_on_vertex(V, C, i, P0)
    
    fi.append(distance.euclidean(force,[0,0,0]))
    
    if distance.euclidean(force,[0,0,0])>5:
        #print(i)
        vi.append(i)
    # if cp.get_perimeter(V, C, i)>4.1:
    #     print(i)
    

# plt.scatter(vi,fi,c="r")


#print(vi)


### cells to look at:
    ## 28,74,76
    
    

cell_list=[28,74,76]


a=cp.get_area(V, C, 3)
p=cp.get_perimeter(V, C, 3)

print(a,p)

pos=nx.get_node_attributes(V,"pos")

length=len(pos)




# x_nodes = [pos[i][0] for i in V]# x-coordinates of nodes
# y_nodes = [pos[i][1] for i in V]# y-coordinates
# z_nodes = [pos[i][2] for i in V]# z-coordinates

# print("set positons")

# x_nodes2 = [pos[i][0] for i in vi]# x-coordinates of nodes
# y_nodes2 = [pos[i][1] for i in vi]# y-coordinates
# z_nodes2 = [pos[i][2] for i in vi]


# edge_list=V.edges()
# x_edges=[]
# y_edges=[]
# z_edges=[]

# #need to fill these with all of the coordiates
# for edge in edge_list:
#     #format: [beginning,ending,None]
#     x_coords = [pos[edge[0]][0],pos[edge[1]][0],None]
    

#     y_coords = [pos[edge[0]][1],pos[edge[1]][1],None]
    

#     z_coords = [pos[edge[0]][2],pos[edge[1]][2],None]
    
#     x_edges += x_coords
#     y_edges += y_coords
#     z_edges += z_coords
 
    
# print("set edges")
# trace_edges = go.Scatter3d(x=x_edges,
#                         y=y_edges,
#                         z=z_edges,
#                         mode='lines',
#                         line=dict(color='black', width=2))

# trace_nodes = go.Scatter3d(x=x_nodes,
#                           y=y_nodes,
#                         z=z_nodes,
#                         mode="markers",
#                         marker=dict(size=2,color="blue"))

# trace_nodes2 = go.Scatter3d(x=x_nodes2,
#                           y=y_nodes2,
#                         z=z_nodes2,
#                         mode="markers",
#                         marker=dict(size=5,color="red"))



# data = [trace_edges, trace_nodes,trace_nodes2]
# fig = go.Figure(data=data)


# fig.show()
# fig.write_html('image.html', auto_open=True)





























## vertices are 109,129, shared cells are 49,50, non shared are 41,57

# f=ff.total_force_on_vertex(V, C, 129, P0)


# ed,ad,pd=ff.energy_deriv_debug(V, C, 109, 41, P0)
#print(ad,pd)



#print(f)

# verts=cp.v_position(V, C, 50)

# for v in verts:
#     plt.scatter(v[0],v[2])



# for i in vxz:
#     plt.scatter(i[0],i[1],c="r")

