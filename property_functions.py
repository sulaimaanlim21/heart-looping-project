# -*- coding: utf-8 -*-
"""
Created on Tue Jan  4 15:05:15 2022

@author: sulai
"""
import numpy as np
import pickle as pkl

import networkx as nx
from scipy.spatial import distance
import random


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

def get_distance(V,C,index1,index2):
    pos=nx.get_node_attributes(V,"pos")
    A,B=pos[index1],pos[index2]
    
    return distance.euclidean(A,B)


def unit_vector(A,B):
    
    Ax,Ay,Az=A[0],A[1],A[2]
    Bx,By,Bz=B[0],B[1],B[2]
    length=distance.euclidean(A,B)
    
    uv=np.array([(Bx-Ax)/length,(By-Ay)/length,(Bz-Az)/length])
    unit_vector=uv
    
    return unit_vector


def get_vector(A,B):
    
    Ax,Ay,Az=A[0],A[1],A[2]
    Bx,By,Bz=B[0],B[1],B[2]
    
    
    
    #length=distance.euclidean([Ax,Ay],[Bx,By])
    
    uv=np.array([(Bx-Ax),(By-Ay),(Bz-Az)])
    vector=uv
    
    return vector


def unit_normal(V,C,A,B,vertex_index,cell_index,n):
    pos=nx.get_node_attributes(V,"pos")
    vertices=nx.get_node_attributes(C,"vertices")
    
    Ax,Ay=A[0],A[1]
    Bx,By=B[0],B[1]
    
    
    
    
    length=distance.euclidean([Ax,Ay],[Bx,By])
    dx=(Ax-Bx)/length
    dy=(Ay-By)/length
   
    
    n1=[-dy,dx]
    n2=[dy,-dx]
    
    ## checking which is outward
    verts=vertices[cell_index]
    
    tot1=0
    tot2=0
    for i in range(len(verts)):
        vert2=verts[i]
        if vert2!=vertex_index:
            Cx1,Cy1=pos[vertex_index][0],pos[vertex_index][1]
            Dx1,Dy1=pos[vert2][0],pos[vert2][1]
            C,D=vector_periodic([Cx1,Cy1], [Dx1,Dy1], n)
            Cx,Cy=C[0],C[1]
            Dx,Dy=D[0],D[1]
            
            
            
            dot1=np.dot(np.array([Dx,Dy])-np.array([Cx,Cy]),n1)
            dot2=np.dot(np.array([Dx,Dy])-np.array([Cx,Cy]),n2)
            tot1+=dot1
            tot2+=dot2
                        
        else:
            pass
    
    
    if tot1<0 and tot2>0:
        n_keep=n1
    elif tot2<0 and tot1>0:
        n_keep=n2
        
    else:
        print("error in outward")
        print(n1,n2)
        n_keep=0
    return n_keep


def get_area(V,C,cell_index,n_size):
 
    
    
   
    ## triangles are 123, 134, 145, 156
    ## A123 = 1/2 mag(12 X 13)
    
    ### form the list of vertices in order:
 
    
    pos=nx.get_node_attributes(V,"pos")
    vertices=nx.get_node_attributes(C,"vertices")
    
    vertex_list=vertices[cell_index]
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
    
        
    total_area=0
    v1=vertex_ordered[0]
    for i in range(len(vertex_ordered)-2):
        j,k=i+1,i+2
        v2,v3=vertex_ordered[j],vertex_ordered[k]
        Cx,Cy=pos[v2][0],pos[v2][1]
        Dx,Dy=pos[v3][0],pos[v3][1]
        
        if distance.euclidean(Cx,pos[v1][0])>n_size/2 and distance.euclidean(Cy,pos[v1][1])<n_size/2:
            if Cx<pos[v1][0]:
                Cx+=n_size*xd
            else:
                Cx-=n_size*xd
                
        elif distance.euclidean(Cx,pos[v1][0])<n_size/2 and distance.euclidean(Cy,pos[v1][1])>n_size/2:
            if Cy<pos[v1][1]:
                Cy+=n_size*yd
            else:
                Cy-=n_size*yd
        elif distance.euclidean(Cx,pos[v1][0])>n_size/2 and distance.euclidean(Cy,pos[v1][1])>n_size/2:
            if Cy<pos[v1][1]:
                Cy+=n_size*yd
            else:
                Cy-=n_size*yd
            if Cx<pos[v1][0]:
                Cx+=n_size*xd
            else:
                Cx-=n_size*xd
    
        else:
            pass
         
        if distance.euclidean(Dx,pos[v1][0])>n_size/2 and distance.euclidean(Dy,pos[v1][1])<n_size/2:
            if Dx<pos[v1][0]:
                Dx+=n_size*xd
            else:
                Dx-=n_size*xd
                
        elif distance.euclidean(Dx,pos[v1][0])<n_size/2 and distance.euclidean(Dy,pos[v1][1])>n_size/2:
            if Dy<pos[v1][1]:
                Dy+=n_size*yd
            else:
                Dy-=n_size*yd
        elif distance.euclidean(Dx,pos[v1][0])>n_size/2 and distance.euclidean(Dy,pos[v1][1])>n_size/2:
            if Dy<pos[v1][1]:
                Dy+=n_size*yd
            else:
                Dy-=n_size*yd
            if Dx<pos[v1][0]:
                Dx+=n_size*xd
            else:
                Dx-=n_size*xd
    
        else:
            pass
     
        area=(1/2)*np.abs(np.cross([[Cx,Cy]-np.array(pos[v1])],[[Dx,Dy]-np.array(pos[v1])]))
        
        total_area+=area

    return total_area


def get_perimeter(V,C,cell_index,n_size):
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
            
            
            
            perimeter+=get_distance(V,C,v,nr,n_size)
            #print(distance.euclidean(v1,v2))
            if get_distance(V,C,v,nr,n_size) >n_size/2:
                print(v,nr)
            else:
                pass
            
            
                    
            
            
        check_list.append(v)
    return perimeter


def area_deriv(V,C,vertex_index,cell_index,n_size):
    pos=nx.get_node_attributes(V,"pos")
    vertices=nx.get_node_attributes(C,"vertices")
    ## edges to look at
    
    neighbours=list(V[vertex_index])
    n_refined=[]
    for i in range(len(neighbours)):
        n_index=neighbours[i]
        if n_index in vertices[cell_index]:
            n_refined.append(n_index)
        else:
            pass
    
    if len(n_refined) !=2:
        print("error in n_refined")
        
    else:
        pass
    
    
    Hijk=vertex_index
    Hgij=n_refined[0]
    Hikl=n_refined[1]
    
    
    
    
    
    Lij=get_distance(V,C,Hijk,Hgij,n_size)
    Lik=get_distance(V,C,Hijk,Hikl,n_size)
    
    
    Nij=unit_normal(V,C,np.array(pos[Hgij]),np.array(pos[Hijk]),vertex_index,cell_index,n_size)
    Nik=unit_normal(V,C,np.array(pos[Hikl]),np.array(pos[Hijk]),vertex_index,cell_index,n_size)
    
    dA_by_dH= (1/2)*(Lij*np.array(Nij)+Lik*np.array(Nik))
    
    return dA_by_dH


def perimeter_deriv(V,C,vertex_index,cell_index,n_size):
    pos=nx.get_node_attributes(V,"pos")
    vertices=nx.get_node_attributes(C,"vertices")
    neighbours=list(V[vertex_index])
    n_refined=[]
    for i in range(len(neighbours)):
        n_index=neighbours[i]
        if n_index in vertices[cell_index]:
            n_refined.append(n_index)
        else:
            pass
    
    if len(n_refined) !=2:
        print("error in n_refined")
    else:
        pass
    
    n_refined_sorted=np.sort(n_refined) ## sort in order for consistency
    
    Hijk=vertex_index
    Hgij=n_refined_sorted[0]
    Hikl=n_refined_sorted[1]
    
    Tij=unit_vector(pos[Hgij],pos[Hijk],n_size)
    Tik=unit_vector(pos[Hikl],pos[Hijk],n_size)
    
    dP_by_dH=np.array(Tij)+np.array(Tik)
    
    return dP_by_dH
            
    
        
def energy_deriv(V,C,vertex_index,cell_index,n_size):
    pos=nx.get_node_attributes(V,"pos")
    vertices=nx.get_node_attributes(C,"vertices")
    area=get_area(V,C,cell_index,n_size)
    perimeter=get_perimeter(V,C,cell_index,n_size)
    dA_by_dH=area_deriv(V,C,vertex_index,cell_index,n_size)
    dP_by_dH=perimeter_deriv(V,C,vertex_index,cell_index,n_size)
    
    dE_by_dH=2*KA*(area-A0)*np.array(dA_by_dH) + 2*KP*(perimeter-P0)*np.array(dP_by_dH)
    
    return dE_by_dH

def energy_deriv_debugging(V,C,vertex_index,cell_index,n_size):
    pos=nx.get_node_attributes(V,"pos")
    vertices=nx.get_node_attributes(C,"vertices")
    area=get_area(V,C,cell_index,n_size)
    perimeter=get_perimeter(V,C,cell_index,n_size)
    dA_by_dH=area_deriv(V,C,vertex_index,cell_index,n_size)
    dP_by_dH=perimeter_deriv(V,C,vertex_index,cell_index,n_size)
    
    dE_by_dH=2*KA*(area-A0)*np.array(dA_by_dH) + 2*KP*(perimeter-P0)*np.array(dP_by_dH)
    
    return dE_by_dH,dA_by_dH,dP_by_dH
        
    
def total_force_on_vertex(V,C,vertex_index,n_size,random_force):
    pos=nx.get_node_attributes(V,"pos")
    vertices=nx.get_node_attributes(C,"vertices")
    cell_list=nx.get_node_attributes(V,"cells")[vertex_index]
    
    forcex=0
    forcey=0
    if random_force==True:
        forcex= random.uniform(0,1)*2
        forcey=random.uniform(0,1)*2
    else:
        pass
        
    for i in range(len(cell_list)):
        cell_index=cell_list[i]
        
        
        forcex-=energy_deriv(V,C,vertex_index,cell_index,n_size)[0]
        forcey-=energy_deriv(V,C,vertex_index,cell_index,n_size)[1]
        
    
    force=np.array([forcex,forcey])
        
   
        
        
    return force

def total_force(V,C,number_of_cells,force_dict,n):
    ##force_dict = {new_list: np.zeros(2,dtype=float) for new_list in V.nodes()}
    pos=nx.get_node_attributes(V,"pos")
    vertices=nx.get_node_attributes(C,"vertices")
    for i in V.nodes():
        
        random_number = random.uniform(0,1)
        if random_number<0.3:
            random_force=True
        else:
            random_force=False
        force=total_force_on_vertex(V,C,i,n,random_force)
        force_dict[i]=np.add(force_dict[i],force)
        

    return force_dict


def d_pos(pos,force,dt,n_size):
    x_new=pos[0]+force[0]*dt
    y_new=pos[1]+force[1]*dt
    if x_new>n_size*xd:
        x_new-=n_size*xd
    elif x_new<0:
        x_new+=n_size*xd
    else:
        pass
    if y_new>n_size*yd:
        y_new-=n_size*yd
    elif y_new<0:
        y_new+=n_size*yd
    else:
        pass
    
    
    return [x_new,y_new]

  
### update all pos
    
def update_all_pos(V,C,force_dict,dt):
    pos=nx.get_node_attributes(V,"pos")
    vertices=nx.get_node_attributes(C,"vertices")
    
    for i in V.nodes():
        V.node[i]["pos"]=d_pos(pos[i],force_dict,dt)
        
        
def position_vertices(V,C,i,vt,size):
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    v_list=[]
    v_ordered=[]
    for v in vertex_ordered:
        v_list.append(pos[v])
    ref_vertice=pos[vt]
    for i in range(len(v_list)):
        test_vertice=v_list[i]
        if distance.euclidean(ref_vertice,test_vertice)>size/2:
            Ax,Ay=ref_vertice[0],ref_vertice[1]
            Bx,By=test_vertice[0],test_vertice[1]
            
            
            if distance.euclidean(Bx,Ax)>size/2 and distance.euclidean(By,Ay)<size/2:
                if Bx<Ax:
                    Bx+=size*xd
                else:
                    Bx-=size*xd
                    
            elif distance.euclidean(Bx,Ax)<size/2 and distance.euclidean(By,Ay)>size/2:
                if By<Ay:
                    By+=size*yd
                else:
                    By-=size*yd
            elif distance.euclidean(Bx,Ax)>size/2 and distance.euclidean(By,Ay)>size/2:
                if By<Ay:
                    By+=size*yd
                else:
                    By-=size*yd
                if Bx<Ax:
                    Bx+=size*xd
                else:
                    Bx-=size*xd
        
            else:
                pass
            vertice_to_add=[Bx,By]
            v_ordered.append(vertice_to_add)
        else:
            v_ordered.append(test_vertice)
            
    return v_ordered