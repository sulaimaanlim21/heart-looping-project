# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 14:26:35 2022

@author: sulai
"""

import numpy as np
import pickle as pkl

import networkx as nx
from scipy.spatial import distance
import random

from tissue_new import tissue_generation
import matplotlib.pyplot as plt 

d=np.sqrt(2/np.sqrt(3))
a=2**(1/2)*3**(-3/4) 
b=2**(-1/2)*3**(-1/4)
A0=1
P0=3.5
KA=1
KP=1
eta=1
xd=d
yd=0.5*3*a





def posi(V,i):
    return V.nodes[i]["pos"]

def get_distance(V,C,index1,index2):
    
    A,B=posi(V,index1),posi(V,index2)
    
    return distance.euclidean(A,B)





def get_vector(A,B):
    
    Ax,Ay,Az=A[0],A[1],A[2]
    Bx,By,Bz=B[0],B[1],B[2]
    
    
    
    #length=distance.euclidean([Ax,Ay],[Bx,By])
    
    uv=np.array([(Bx-Ax),(By-Ay),(Bz-Az)])
    vector=uv
    
    return vector


def unit_vector(V,a,b):
    A,B=posi(V,a),posi(V,b)
    
    length=distance.euclidean(A,B)
    uv=np.array(B)/length-np.array(A)/length
    return uv
    
    
def get_perimeter(V,C,cell_index):
    pos=nx.get_node_attributes(V,"pos")
    vertices=nx.get_node_attributes(C,"vertices")
    vertex_list=vertices[cell_index]
    perimeter=0
    check_list=[]
    # print(len(vertex_list))
    # print(vertex_list)
    checking=False
    
    for i in range(len(vertex_list)):
        v=vertex_list[i]
        
        neighbours=list(V[v])
        n_refined=[]
        
        for j in range(len(neighbours)):
            n=neighbours[j]
            if n in check_list:
                pass
            
            
            
            elif n in vertex_list:
                
                n_refined.append(n)
            else:
                pass
            
        
        for k in range(len(n_refined)):
            nr=n_refined[k]
            
            
            
            perimeter+=distance.euclidean(pos[v],pos[nr])
            #print(distance.euclidean(v1,v2))
            
            
            
        check_list.append(v)
    return perimeter

def order(V,C,i):
    vertice_indices=C.nodes[i]["vertices"]
    pos=nx.get_node_attributes(V,"pos")
    
    
    vertices=nx.get_node_attributes(C,"vertices")
    
    vertex_list=vertice_indices
    vertex_ordered=[vertex_list[0]]
    neighbours0=list(V[vertex_list[0]])
    n_refined0=[]
    for neighbour0 in neighbours0:
        if neighbour0 in vertex_list:
            
            n_refined0.append(neighbour0)
        else:
            pass
        
    neighbour1=n_refined0[0]
    vertex_ordered.append(neighbour1)
    ### set up first 2, now use iteration to set up rest
    
    while len(vertex_ordered)<len(vertex_list):
        vert=vertex_ordered[-1]
        neighbours=list(V[vert])
        
        for neighbour in neighbours:
            if neighbour in vertex_list:
                if neighbour not in vertex_ordered:
                    vertex_ordered.append(neighbour)
                else:
                    pass
                
            else:
                pass
    
    

    
    return vertex_ordered

def order_rh(V,C,i):
    vertice_indices=C.nodes[i]["vertices"]
    pos=nx.get_node_attributes(V,"pos")
    
    
    vertices=nx.get_node_attributes(C,"vertices")
    
    vertex_list=vertice_indices
    vertex_ordered=[vertex_list[0]]
    neighbours0=list(V[vertex_list[0]])
    n_refined0=[]
    for neighbour0 in neighbours0:
        if neighbour0 in vertex_list:
            
            n_refined0.append(neighbour0)
        else:
            pass
        
    neighbour1=n_refined0[0]
    vertex_ordered.append(neighbour1)
    ### set up first 2, now use iteration to set up rest
    
    while len(vertex_ordered)<len(vertex_list):
        vert=vertex_ordered[-1]
        neighbours=list(V[vert])
        
        for neighbour in neighbours:
            if neighbour in vertex_list:
                if neighbour not in vertex_ordered:
                    vertex_ordered.append(neighbour)
                else:
                    pass
                
            else:
                pass
    
    
    ##check to see if right handed (unit normal points outwards)
    va,vb=posi(V,vertex_ordered[0]),posi(V,vertex_ordered[1])
    centroid=get_centroid(V,C,i)
    
    v1,v2=va-centroid,vb-centroid
    v_n=np.cross(v1,v2)
    ## check if outwards:
    v_ref=C.nodes[i]["normal"]
    dot=np.dot(v_n,v_ref)
    if dot>0:
        pass
    else:
        vertex_ordered.reverse()
        
    
    return vertex_ordered





def v_position(V,C,i):
    index_ordered=order(V,C,i)
    v_ordered=[posi(V,j) for j in index_ordered]
    
    return v_ordered




def get_area(V,C,cell_index):
 
    
    
   
    ## triangles are 123, 134, 145, 156
    ## A123 = 1/2 mag(12 X 13)
    
    ### form the list of vertices in order:
        
    vertices=v_position(V,C,cell_index)
    
    
    
    
    
    
    v_centroid=get_centroid(V, C, cell_index)
    
    #print("vertices are", vertices)
    #print("centroid is", v_centroid)
    nv=len(vertices)
    total_area=0
    for i in range(nv-1):
        ##triangle a,b,c where c is centroid
        a,b=vertices[i],vertices[i+1]
        
        e1=np.array(a)-np.array(v_centroid)
        e2=np.array(b)-np.array(v_centroid)
        #print("e1,e2 are", e1,e2)
        
        #print(np.cross(e1,e2),np.abs(np.cross(e1,e2)))
        area=(1/2)*distance.euclidean((np.cross(e1,e2)),[0,0,0])
        #print("area of tri is", area)
        total_area+=area
                
    ## vert[-1] with vert[0]   
    
    a,b=vertices[0],vertices[-1]
    e1=np.array(a)-np.array(v_centroid)
    e2=np.array(b)-np.array(v_centroid)

    area=(1/2)*distance.euclidean((np.cross(e1,e2)),[0,0,0])
    total_area+=area
    
    
    return total_area



def unit_normal(V,C,a,b,cell_index):
    ## cell index is i
    ## vectors in order cabd 
    
    ##get c
    
    vert_indices=C.nodes[cell_index]["vertices"]
    
    for j in list(V[a]):
        if j in vert_indices and j!=b:
            c=j
        else:
            pass
        
    for k in list(V[b]):
        if k in vert_indices and k!=a:
            d=k
        else:
            pass
        
    ## plane vectors 1 and 2, pv1 and pv2:
        
    va,vb=posi(V,a),posi(V,b)
    vc,vd=posi(V,c),posi(V,d)
    
    pv1=np.array(va)-np.array(vc)
    pv2=np.array(vb)-np.array(vd)
    
    
    ## plane normal, pn:
    pn=np.cross(pv1,pv2)
    
    
    vec=np.array(va)-np.array(vb)
    
    
    ##unit normal is cross of normal and vec:
    
    normal=np.cross(pn,vec)
    
    un_length=distance.euclidean(normal,[0,0,0])
    un=(1/un_length)*np.array(normal)
    
    
    ##check if points outwards
    
    
    v_centroid=get_centroid(V, C, cell_index)
    
    
    midpoint = 0.5*np.array(va)+0.5*np.array(vb)
    mid_centroid=np.array(midpoint)-np.array(v_centroid)
    
    dot_product=np.dot(mid_centroid,un)
    
    if dot_product>0:
        unit_normal=np.array(un)
    else:
        unit_normal=-np.array(un)
        
        
        
    return unit_normal
    
 
def get_centroid(V,C,cell_index):
    vertices=v_position(V,C,cell_index)
    
    centroid=np.array([0,0,0])
    length=len(vertices)
    for vert in vertices:
        centroid=centroid+(1/length)*np.array(vert)
    
    return centroid
 
    
def cell_normal_initial(V,C,cell_index):
    centroid=get_centroid(V,C,cell_index)
    vert_list=v_position(V,C,cell_index) ## vertices positions in order
    v_tot=np.array([float(0),float(0),float(0)])
    for i in range(len(vert_list)):
        a=vert_list[i]
        b=vert_list[(i+1)%len(vert_list)]
        v_a,v_b=np.array(a)-np.array(centroid),np.array(b)-np.array(centroid)
        
        v_n=np.cross(v_a,v_b)
        mag=distance.euclidean(v_n,[0,0,0])
        v_un=(1/mag)*np.array(v_n)
        
        ## check if faces outward
        v_check=np.array([v_un[0],v_un[1]])
        v_2check=np.array([centroid[0],centroid[1]])
        if np.dot(v_check,v_2check)>0:
            v=v_un
        else:
            v=-v_un
            
        v_tot+=v
        
    
    v_tot_mag=distance.euclidean(v_tot,[0,0,0])
    v_tot_un=(1/v_tot_mag)*v_tot
    
    
    return v_tot_un
            
        
        
        
        
        
        
        
def cell_normal(V,C,cell_index):
    centroid=get_centroid(V,C,cell_index)
    vert_list=v_position(V,C,cell_index) ## vertices positions in order
    v_tot=np.array([float(0),float(0),float(0)])
    for i in range(len(vert_list)):
        a=vert_list[i]
        b=vert_list[(i+1)%len(vert_list)]
        v_a,v_b=np.array(a)-np.array(centroid),np.array(b)-np.array(centroid)
        
        v_n=np.cross(v_a,v_b)
        mag=distance.euclidean(v_n,[0,0,0])
        v_un=(1/mag)*np.array(v_n)
        
        ## check if faces outward
        v_check=np.array(v_un)
        v_check2=C.nodes[cell_index]["normal"]
        if np.dot(v_check,v_check2)>0:
            v=v_un
        else:
            v=-v_un
            
        v_tot=v_tot+np.array(v)
        
    
    v_tot_mag=distance.euclidean(v_tot,[0,0,0])
    v_tot_un=(1/v_tot_mag)*v_tot
    
    
    return v_tot_un
  
def set_initial_cell_normals(V,C):
    for c in C:
        C.nodes[c]["normal"]=cell_normal_initial(V,C,c)
        
def set_cell_normals(V,C):
    for c in C:
        C.nodes[c]["normal"]=area_properties(V, C, c)[2]
        
        
def area_perimeter_calculate(V,C):
    for c in C:
        area,area_vec,unit_normal=area_properties(V, C, c)
        C.nodes[c]["area"]=area
        C.nodes[c]["area_vector"]=area_vec
        C.nodes[c]["perimeter"]=get_perimeter(V,C,c)
        
def set_centroids(V,C):
    for c in C:
        C.nodes[c]["centroid"]=get_centroid(V, C, c)
        
def set_cell_properties(V,C):
    set_cell_normals(V, C)
    area_perimeter_calculate(V,C)
    set_centroids(V, C)
    
def set_initial_cell_properties(V,C):
    set_initial_cell_normals(V, C)
    area_perimeter_calculate(V,C)
    set_centroids(V, C)
    
    
    
    
def area_properties(V,C,c):
    vertex_list=order_rh(V, C, c)
    
    area_total=np.array([0,0,0])
    n=len(vertex_list)
    
    for i in range(n):
        
        a_i=i
        b_i=(i+1) % n
        
        
    
        
        a=vertex_list[a_i]
        b=vertex_list[b_i]
        v_area=0.5*np.cross(posi(V,a),posi(V,b))
        
        area_total=area_total+v_area
        
    
    area_mag=distance.euclidean([0,0,0],area_total)
    area_vec=area_total
    
    unit_normal=(1/area_mag)*area_vec
    return area_mag,area_vec,unit_normal




def area_deriv_total(V,C,v):
    cell_list=V.nodes[v]["cells"]
    area_deriv_total=np.array([0,0,0])
    for c in cell_list:
        area_deriv=area_deriv_c(V, C, v, c)
        area_deriv_total=area_deriv_total+area_deriv
        
    return area_deriv_total
        

def area_deriv_c(V,C,v,c):
    ## area deriv is (1/A)*sum(Ak *dAk/dx)
    ## dAk/dx = (1/2)*(x_[i-1]-x_[i+1]) X ek    (ek is unit vector in k)
    ##call x+1 index b x-1 index a
    
    ex,ey,ez=np.array([1,0,0]),np.array([0,1,0]),np.array([0,0,1])
    A_mag,A_vec=C.nodes[c]["area"],C.nodes[c]["area_vector"]
    vertex_list=order_rh(V, C, c)
    n=len(vertex_list)
    i=vertex_list.index(v)
    b=vertex_list[(i+1)%n]
    a=vertex_list[i-1]
    
    
    dAx=(1/2)*A_vec[0]*(np.cross(np.array(posi(V,b))-np.array(posi(V,a)),ex))
    dAy=(1/2)*A_vec[1]*(np.cross(np.array(posi(V,b))-np.array(posi(V,a)),ey))
    dAz=(1/2)*A_vec[2]*(np.cross(np.array(posi(V,b))-np.array(posi(V,a)),ez))
    
    dA_by_dx=(1/A_mag)*(dAx+dAy+dAz)
    
    return np.array(dA_by_dx)
