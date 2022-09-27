# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 14:43:02 2022

@author: sulai
"""

import numpy as np
import pickle as pkl
import networkx as nx
from scipy.spatial import distance
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
import cell_properties as cp
import force_functions_new as ff
pio.kaleido.scope.mathjax = None
step_number=0

k=500
path='C:/Users/sulai/OneDrive - Imperial College London/PhD/Simulation/files/'
file_name = 'stage_new2'
file_name2 = 'stagec_new2'
V = nx.read_gpickle(path+file_name + '.pickle')
pos=nx.get_node_attributes(V,"pos")
C=nx.read_gpickle(path+file_name2 + '.pickle')

# v_p=pos[404]
# V.nodes[404]["pos"]=[v_p[0]+0.5,v_p[1],v_p[2]]

# for v in V:
#     if distance.euclidean(cp.posi(V,v),[2.7,-1.1,13.4])<0.1:
#         print(v)




pos=nx.get_node_attributes(V,"pos")
## node is 503

bending_force=ff.bending_force_new(V,C,503)
mag=distance.euclidean(bending_force,[0,0,0])
unit_vector=(1/mag)*np.array(bending_force)

pos_bf=pos[503]+unit_vector




















length=len(pos)
for i in V:
    V.nodes[i]["checked"]=False


        
        
        
        
           




x_nodes = [pos[i][0] for i in V ]# x-coordinates of nodes
y_nodes = [pos[i][1] for i in V ]# y-coordinates
z_nodes = [pos[i][2] for i in V ]# z-coordinates

x_nodes2 = [pos_bf[0]] # x-coordinates of nodes
y_nodes2 = [pos_bf[1]] # y-coordinates
z_nodes2 = [pos_bf[2]]# z-coordinates



# x_nodes2 = [V.nodes[i]["anchor_pos"][0] for i in V if V.nodes[i]["anchor"]==True]# x-coordinates of nodes
# y_nodes2 = [V.nodes[i]["anchor_pos"][1] for i in V if V.nodes[i]["anchor"]==True]# y-coordinates
# z_nodes2 = [V.nodes[i]["anchor_pos"][2] for i in V if V.nodes[i]["anchor"]==True]# z-coordinates


  
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
    
    
    
    # x_edges2 += x_coords
    # y_edges2 += y_coords
    # z_edges2 += z_coords

    x_edges += x_coords
    y_edges += y_coords
    z_edges += z_coords
        
 
    

trace_edges = go.Scatter3d(x=x_edges,
                        y=y_edges,
                        z=z_edges,
                        mode='lines',
                        line=dict(color='black', width=2))

# trace_edges2 = go.Scatter3d(x=x_edges2,
#                         y=y_edges2,
#                         z=z_edges2,
#                         mode='lines',
#                         line=dict(color='blue', width=1))

trace_nodes = go.Scatter3d(x=x_nodes,
                          y=y_nodes,
                        z=z_nodes,
                        mode="markers",
                        marker=dict(size=3,color="red")
                        )

trace_nodes2 = go.Scatter3d(x=x_nodes2,
                          y=y_nodes2,
                        z=z_nodes2,
                        mode="markers",
                        marker=dict(size=3,color="red")
                        )



data = [trace_edges,trace_nodes2]

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
# layout = Layout(
#     paper_bgcolor='rgba(0,0,0,0)',
#     plot_bgcolor='rgba(0,0,0,0)'
# )

# fig['layout']['paper_bgcolor']='rgba(0,0,0,0)'
# fig['layout']['plot_bgcolor']='rgba(0,0,0,0)'

fig['layout']['scene']['aspectmode']='data'
fig.update_layout(scene_camera=camera)
#fig.update_layout(
# scene = dict(
#     xaxis = dict(range=[-8,8],),
#                  yaxis = dict( range=[-8,8],),
#                  zaxis = dict( range=[-2,20],)))

# fig.update_xaxes(visible=False) # Hide x axis ticks 
# fig.update_yaxes(visible=False)
# #fig.update_zaxes(visible=False)
# # data = [trace_edges, trace_nodes]
# # fig = go.Figure(data=data)
# fig.update_layout(title={
#                     'text': "Hexagonally tiled 3D tube",
#                     'y':1,
#                     'x':0.5,
#                     },
#                   scene = dict(
#                     xaxis = dict(
#                          backgroundcolor="rgb(0, 0, 0)",
#                          gridcolor="white",
#                          showbackground=False,
#                          zerolinecolor="white",
#                          visible=False),
#                     yaxis = dict(
#                         backgroundcolor="rgb(0, 0, 0)",
#                          gridcolor="white",
#                          showbackground=False,
#                          zerolinecolor="white",
#                          visible=False),
#                     zaxis = dict(
#                         backgroundcolor="rgb(0, 0, 0)",
#                          gridcolor="white",
#                          showbackground=False,
#                          zerolinecolor="white",
#                          visible=False),)
                   
#                   )



#fig.show()
fig.write_html('image.html', auto_open=True)





#pio.write_image(fig, 't' + str(step_number)+'.png', format=None,
                  #scale=None, width=None, height=None)


#fig.write_image('t' + str(step_number)+'.png',engine='orca')
#plt.savefig(file_name + '.png')

#plt.clf()
      
    