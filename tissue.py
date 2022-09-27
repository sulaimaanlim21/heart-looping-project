# -*- coding: utf-8 -*-
"""
Created on Wed Dec 15 16:41:50 2021

@author: sulai
"""


import numpy as np
import pickle as pkl
from cell_creation import generate_cells
import networkx as nx
from scipy.spatial import distance
from hexalattice.hexalattice import *
from mpl_toolkits.mplot3d import Axes3D
#from property_functions import total_force_on_vertex
d=np.sqrt(2/np.sqrt(3))
a=2**(1/2)*3**(-3/4) ## distance from center to edge and side length
b=2**(-1/2)*3**(-1/4) ## center to flat edge distance



def tissue_generation(n_x,n_y):
    cells=generate_cells(n_x,n_y)
    
    x_displacement=n_x*d
    xd=x_displacement
    yd=(n_y/2)*3*a
    #print("cells initialised")
    
    ##get max and min locations
    
    
    
    vertices_list=[]
    for i in range(len(cells)):
        vertices=cells[i]["vertices"]
        for j in range(6):
            vertices_list.append(vertices[j])
            
    ### reduce list
    
    vrefined=vertices_list[0:6]
    
    v_rest=vertices_list
    
    for i in range(len(v_rest)):
       add=True
       for j in range(len(vrefined)):
           
           if distance.euclidean(vrefined[j],v_rest[i]) < a/2:
               add=False
               break
       
       if add==True:
           vrefined.append(v_rest[i])
       else:
           pass
            
    
    
    ### anchor top and bottom vertices 
    
    
    
    
    
    V=nx.Graph()
    
    for i in range(len(vrefined)):
        V.add_node(i,pos=vrefined[i])
        
        
        if distance.euclidean(yd,vrefined[i][1])<3*a/4:
            V.nodes[i]["anchor"]=True
        elif distance.euclidean(0,vrefined[i][1])<a:
            V.nodes[i]["anchor"]=True
        else:
            V.nodes[i]["anchor"]=False
            
       
        
    
     
        
        
    ## Graph generated for vertices
        

    
        
    pos=nx.get_node_attributes(V, "pos")    
        
    for i in V:
        for j in V:
            if distance.euclidean(pos[i],pos[j])<a+0.1 and i!=j:
                V.add_edge(i,j)
            else:
                pass
            
            
    for i in V:
        cell_list=[]
        for j in range(len(cells)):
            center=cells[j]["center"]
            if distance.euclidean(pos[i],center)<a+0.1:
                cell_list.append(j)
            else:
                pass
        V.nodes[i]["cells"]=cell_list
        
    
    list_to_remove_x=[]
    for i in V:
        posx_shift,posy=pos[i][0]+xd,pos[i][1]
        for j in V:
            if distance.euclidean([posx_shift,posy],pos[j])<0.01:
                neighbours=list(V[j])
                cells=V.nodes[j]["cells"]
                for k in range(len(neighbours)):
                    V.add_edge(i,neighbours[k])
                for l in range(len(cells)):
                    
                    cell_list=V.nodes[i]["cells"]
                    cell_list.append(cells[l])
                    
                    V.nodes[i]["cells"]=cell_list
                    
                
                
                    
                list_to_remove_x.append(j)
            else:
                pass
 
    
    
    
            
    for i in range(len(list_to_remove_x)):
     
      index=list_to_remove_x[i]
      if index in V:
         
          V.remove_node(index)
      else:
          pass
            
            
            
    cells=generate_cells(n_x,n_y)
    C=nx.Graph()

    for i in range(len(cells)):
        v_list=[]
        for j in V:
            if i in V.nodes[j]["cells"]:
                v_list.append(j)
            else:
                pass
      
        
                
        C.add_node(i,vertices=v_list)  
        
        
    
    for i in C:
        if len(C.nodes[i]["vertices"])!=6:
            
            v_list=C.nodes[i]["vertices"]
            for i in range(len(v_list)):
                plt.scatter(pos[v_list[i]][0],pos[v_list[i]][1],c="r")
            
            
        else:
            pass
        
    # for i in V:
    #     if len(V.nodes[i]["cells"])!=3:
    #         plt.scatter(pos[i][0],pos[i][1],c="g")
    #     else:
    #         pass
        
    # check2=0
    # for i in V:
    #     neighbours=list(V[i])
    #     if len(neighbours)!=3:
    #         check2+=1
    #         print("vertex that does not have 3 neighbours is", i, "number of neighbours is", len(neighbours))
    #     else:
    #         pass
    # if check2!=0:
    #     print("error, not all vertices have 3 neighbours")
    
    
    ### check for associated cells
    
    # check1=0
    # for i in V:
        
    #     cells_list=V.nodes[i]["cells"]
    #     if len(cells_list)!=3:
    #         check1+=1
    #         print("vertex", i, "has", len(cells_list), "associated cells which are", cells_list)
    #     else:
    #         pass

    
    
        
    # for i in C:
    #     if len(C.nodes[i]["vertices"])!=6:
            
    #         print("error, cell", i, "has", len(C.nodes[i]["vertices"]), "vertices" )
            
    #     else:
    #         pass 
    
    
    
    ## making 3D:
    for i in V:
        
        ##r = C/2*pi
        r=xd/(2*np.pi)
        
        x_c = r*np.cos(2*np.pi*pos[i][0]/xd) # x-coordinates of nodes
        y_c = r*np.sin(2*np.pi*pos[i][0]/xd)# y-coordinates
        z_c = pos[i][1] # z-coordinates
        V.nodes[i]["pos"]=[x_c,y_c,z_c]



    ## making 3D:
    
        
    mid=n_y*a*(3/4)
    lb=mid-1*a
    ub=mid+2*a

    ## add in the shrink

    

    
                
    
    
    return V,C


# V,C=tissue_generation(12,16)
# #print(len(V))
# pos=nx.get_node_attributes(V,"pos")
# plot=nx.draw_networkx(V,pos,node_size=1,node_color="blue",node_shape="+",with_labels=False)

# for i in V:
#     if V.nodes[i]["anchor"]==True:
#         plt.scatter(pos[i][0],pos[i][1],c="r",marker="x")



# # # 3d spring layout

# # Extract node and edge positions from the layout
# node_xyz = np.array([pos[v] for v in sorted(V)])
# edge_xyz = np.array([(pos[u], pos[v]) for u, v in V.edges()])

# # Create the 3D figure
# fig = plt.figure()
# ax = fig.add_subplot(111, projection="3d")

# # Plot the nodes - alpha is scaled by "depth" automatically
# ax.scatter(*node_xyz.T, s=10, ec="w",marker="x")

# # Plot the edges
# for vizedge in edge_xyz:
#     ax.plot(*vizedge.T, color="tab:gray")


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


# _format_axes(ax)
# fig.tight_layout()
# plt.show()
