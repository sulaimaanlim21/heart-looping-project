# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 14:09:30 2022

@author: sulai
"""




import numpy as np



import networkx as nx
from scipy.spatial import distance
from cell_properties import v_position,get_centroid
#from shapely.geometry import Polygon



### if edge between 1 an 3 is less than Lt, transition occurs, edge forms between 2 and 4
### need position of vertices between 1 and 3, move them and form new edges. Need centroid
##start with vectors v1 and v2

def t1_transition(V,C,a,b,lt):
    pos=nx.get_node_attributes(V,"pos")
    
    
    #lt=0.1
    #print("a cells is", V.nodes[a]["cells"])
    cell_list=[]
    for i in range(len(V.nodes[a]["cells"])):
        cell_list.append(V.nodes[a]["cells"][i])
        
    
    for i in range(len(V.nodes[b]["cells"])):
        cell=V.nodes[b]["cells"][i]
        if cell not in cell_list:
            cell_list.append(cell)
        else:
            pass
    
       
    c24=[]
    
    for i in range(len(cell_list)):
        cell=cell_list[i]
        if cell not in V.nodes[a]["cells"]:
            c3=cell
        elif cell not in V.nodes[b]["cells"]:
            c1=cell
        elif cell in V.nodes[a]["cells"] and cell in V.nodes[b]["cells"]:
            c24.append(cell)
        else:
            print("error in asigning cells for t1")
    #print("cell list is", cell_list)
    
    c2=c24[0]
    c4=c24[1]
    
    ##Identify centre of the system: ### consider periodicity!!!
   
            
            
    
    midpoint = 0.5*(np.array(pos[a])+np.array(pos[b]))
    
    
    ##centroids
    
    
    
    centroid2=get_centroid(V, C, c2)
    
    
    
    centroid4=get_centroid(V,C,c4)
    
    
    
    
    centroid_vector=np.array(centroid4)-np.array(centroid2)
    mag=distance.euclidean(centroid_vector,[0,0,0])
    
    centroid_unit_vector=(1/mag)*centroid_vector
    print(distance.euclidean(centroid_unit_vector,[0,0,0]))
    position_a=midpoint+0.5*1.5*lt*centroid_unit_vector
    position_b=midpoint-0.5*1.5*lt*centroid_unit_vector
    
    ##update position:
    #print("previous pos was ",V.nodes[a]["pos"], "new pos is", position_a)
    ##make periodic
    
       
    V.nodes[a]["pos"]=position_a
    V.nodes[b]["pos"]=position_b
    
    ### update edges and cell neighbours:
    
    ##identify alpha and delta:
    a_neighbours=list(V[a])
    for neighbour in a_neighbours:
        
        if neighbour in C.nodes[c1]["vertices"] and neighbour in C.nodes[c4]["vertices"]:
            v_alpha=neighbour
        elif neighbour in C.nodes[c1]["vertices"] and neighbour in C.nodes[c2]["vertices"]:
            v_delta=neighbour
        else:
            pass 
        
    b_neighbours=list(V[b])
    #print(b_neighbours,C.nodes[c3]["vertices"],C.nodes[c4]["vertices"])
    for neighbour in b_neighbours:
        
        if neighbour in C.nodes[c3]["vertices"] and neighbour in C.nodes[c4]["vertices"]:
            v_beta=neighbour
        elif neighbour in C.nodes[c2]["vertices"] and neighbour in C.nodes[c3]["vertices"]:
            v_gamma=neighbour
        else:
            pass 
                
    ## remove and replace edges, replace cell list for vertices, replace vertice list for cells
    #print(v_alpha,v_delta)
   # print(a_neighbours,C.nodes[c1]["vertices"],C.nodes[c4]["vertices"])
    ##edges:
    V.remove_edge(a,v_delta)
    V.add_edge(a,v_beta)
    
    V.remove_edge(b,v_beta)
    V.add_edge(b,v_delta)
    
    ##cell list
    cell_list_a=[c1,c3,c4]
    cell_list_b=[c1,c2,c3]
    
    V.nodes[a]["cells"]=cell_list_a
    V.nodes[b]["cells"]=cell_list_b
    
    ##vertices_lists
    
    c1v_list=C.nodes[c1]["vertices"]
    c2v_list=C.nodes[c2]["vertices"]
    c3v_list=C.nodes[c3]["vertices"]
    c4v_list=C.nodes[c4]["vertices"]
    
    c1v_list.append(b)
    c3v_list.append(a)
    
    c2v_list.remove(a)
    c4v_list.remove(b)
    # 
    
