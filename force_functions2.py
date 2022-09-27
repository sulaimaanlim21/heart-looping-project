# -*- coding: utf-8 -*-
"""
Created on Tue Feb  1 13:43:07 2022

@author: sulai
"""

import numpy as np
import pickle as pkl

import networkx as nx
from scipy.spatial import distance
import random
from shapely.geometry import Polygon
from tissue import tissue_generation
import matplotlib.pyplot as plt 

import cell_properties as cp
from cell_properties import posi

KA=1
KP=1
eta=1
dt=0.01
#lt=0.015
A0=1
A1=1
P1=3.4*np.sqrt(A1)
a=2**(1/2)*3**(-3/4) 
P0=3.4
A2=1.4
P2=3.4*np.sqrt(A2)
A3=1.6
P3=3.4*np.sqrt(A3)

mid=12*a*(3/4)

def energy_deriv(V,C,vertex_index,cell_index,parameters):
    A1,A2,A3=parameters[0],parameters[1],parameters[2]
    if C.nodes[cell_index]["type"]=="avc":
        gm=1
    else:
        gm=growth_factor3(V,C,cell_index)
    # if C.nodes[cell_index]["region"]=="ul":
    #     gm=1.5
    # elif C.nodes[cell_index]["region"]=="ur":
    #     gm=1
    # elif C.nodes[cell_index]["region"]=="ll":
    #     gm=1
    # else:
        
    #     gm=1.5
        
    if C.nodes[cell_index]["type"]=="atr":
        A,P=A3*gm,3.4*np.sqrt(A3*gm)
        
    elif C.nodes[cell_index]["type"]=="avc":
        angle=C.nodes[cell_index]["angle"] 
        
        A,P= A1*gm,3.4*np.sqrt(A1*gm)
    else:
        A,P= A2*gm,3.4*np.sqrt(A2*gm)
    area=cp.get_area(V,C,cell_index)
    perimeter=cp.get_perimeter(V,C,cell_index)
    dA_by_dH=area_deriv(V,C,vertex_index,cell_index)
    dP_by_dH=perimeter_deriv(V,C,vertex_index,cell_index)
    
    dE_by_dH=2*(area-A)*np.array(dA_by_dH) + 2*(perimeter-P)*np.array(dP_by_dH)
    
    return dE_by_dH


def grow_angle(V,C,i,mid):
    
    
    p=cp.get_centroid(V, C, i)
    pz=p[2]
    yd=distance.euclidean(mid,pz)
      
          
    x_array=[]
    y_array=[]
    for cell in C:
        p1=cp.get_centroid(V, C, cell)
        if pz-a/2<p1[2]<pz+a/2:
            x_array.append(p1[0])
            y_array.append(p1[1])
        else:
            pass
                  
    xc,yc=np.mean(x_array),np.mean(y_array)
             
    
    ##get angle
    x,y,z=p[0]-xc,p[1]-yc,p[2]
    radius=np.sqrt(x**2+y**2)
    ### set angle correctly
    if y>0:
        angle=np.arccos(x/radius)
    else:
        angle=2*np.pi-np.arccos(x/radius)
        
    C.nodes[i]["angle"]=angle
    #print('set angle as',angle)
    
def grow_factor_z(V,C,i,mid):
    p=cp.get_centroid(V, C, i)
    pz=p[2]
    yd=distance.euclidean(mid,pz)
    growth_factor=0.5*np.tanh(2*np.pi*(pz-mid)/mid)
    C.nodes[i]["growth_z"]=growth_factor      
    
    
def growth_factor(V,C,i):
    growth_z=C.nodes[i]["growth_z"]
    #if cp.get_centroid(V,C,i)[2]<mid:
    ref_angle=np.pi*(1/2)
    # else:
    #     ref_angle=np.pi*(1/2)
    angle=C.nodes[i]["angle"]   
    if np.pi>angle>0:
        
        growth_factor=0.5+1/(np.cosh(ref_angle-angle))#-growth_z
    else:
        growth_factor=1/(np.cosh(ref_angle-angle))#+growth_z
    
    # if growth_factor>0:
    growth_modifier=1+growth_factor
    # else:
    #     growth_modifier=1+0.2*growth_factor
    
    return growth_modifier

def growth_factor3(V,C,i):
    growth_z=C.nodes[i]["growth_z"]
    if cp.get_centroid(V,C,i)[2]<mid:
        ref_angle=np.pi*(1/2)
    else:
        ref_angle=np.pi*(3/2)
    angle=C.nodes[i]["angle"]   
    
        
    growth_factor=1/(2*(np.cosh(ref_angle-angle)))#-growth_z
    
    # if growth_factor>0:
    growth_modifier=1+growth_factor
    # else:
    #     growth_modifier=1+0.2*growth_factor
    
    return growth_modifier






def area_deriv(V,C,vertex_index,cell_index):
    
    
    ## edges to look at
    
    neighbours=list(V[vertex_index])
    n_refined=[]
    for n_index in neighbours:
        
        if n_index in C.nodes[cell_index]["vertices"]:
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
    
    
    
    
    
    Lij=cp.get_distance(V,C,Hijk,Hgij)
    Lik=cp.get_distance(V,C,Hijk,Hikl)
    
    
    Nij=cp.unit_normal(V,C,Hgij,Hijk,cell_index)
    Nik=cp.unit_normal(V,C,Hikl,Hijk,cell_index)
    
    dA_by_dH= (1/2)*(Lij*np.array(Nij)+Lik*np.array(Nik))
    
    return dA_by_dH


def perimeter_deriv(V,C,vertex_index,cell_index):
    
    neighbours=list(V[vertex_index])
    n_refined=[]
    for n_index in neighbours:
        
        if n_index in C.nodes[cell_index]["vertices"]:
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
    
    Tij=cp.unit_vector(V,Hgij,Hijk)
    Tik=cp.unit_vector(V,Hikl,Hijk)
    
    dP_by_dH=np.array(Tij)+np.array(Tik)
    
    return dP_by_dH




def energy_deriv_debug(V,C,vertex_index,cell_index):
    # if C.nodes[cell_index]["type"]=="atr":
    #     A,P= A0,P0
    # elif C.nodes[cell_index]["type"]=="avc":
    #     A,P= A1,P1
    # else:
    A,P= A3,P3
    area=cp.get_area(V,C,cell_index)
    perimeter=cp.get_perimeter(V,C,cell_index)
    dA_by_dH=area_deriv(V,C,vertex_index,cell_index)
    dP_by_dH=perimeter_deriv(V,C,vertex_index,cell_index)
    
    dE_by_dH=2*(area-A)*np.array(dA_by_dH) + 2*(perimeter-P)*np.array(dP_by_dH)
    
    return dE_by_dH,dA_by_dH,dP_by_dH,area,perimeter



def total_force_on_vertex(V,C,vertex_index,parameters):
    
    cell_list=V.nodes[vertex_index]["cells"]
    
    
    total_force=np.array([0,0,0])
        
    for cell_index in cell_list:
        
        
        force=-energy_deriv(V,C,vertex_index,cell_index,parameters)
        total_force=total_force+np.array(force)
        
    return total_force



    

def pressure_force(V,v,r):
    x,y=cp.posi(V,v)[0],cp.posi(V,v)[1]
    
    ## make radial unit vector
    mag=distance.euclidean([x,y],[0,0])
    force_vec=(1/mag)*np.array([x,y,0])
    
    
    ## force magnitude
    force_mag=np.exp(r-mag)
    
    pressure_force=force_mag*force_vec
    
    return pressure_force


def pressure_force2(V,i,r,mag):
    
    
        
    pz=posi(V,i)[2]
    
    
    
    x_array=[]
    y_array=[]
    for p1 in V:
        if pz-a/2<posi(V,p1)[2]<pz+a/2:
            x_array.append(posi(V,p1)[0])
            y_array.append(posi(V,p1)[0])
        else:
            pass
                    
    xc,yc=np.mean(x_array),np.mean(y_array)
    x,y,z=posi(V,i)[0]-xc,posi(V,i)[1]-yc,posi(V,i)[2]
    radius=np.sqrt(x**2+y**2)
    
    ### set angle correctly
    if y>0:
        angle=np.arccos(x/radius)
    else:
        angle=2*np.pi-np.arccos(x/radius)
        
    
    
    
    
    
    
    
    
    ## make radial unit vector
    
    force_vec=np.array([np.cos(angle),np.sin(angle),0])
    
    
    ## force magnitude
    force_mag=mag*np.exp(r-radius)
    
    pressure_force=force_mag*force_vec
    
    return pressure_force


def bending_force(V,C,i,k):
    #### need the centre of an edge, and the cell centre
    
    tot_force=np.array([0,0,0])
    if V.nodes[i]["anchor"]==False:
        for n in list(V[i]):
            ## define C1,C2,r
            r= 0.5*np.array(cp.posi(V,i))+0.5*np.array(cp.posi(V,n))
            cells=[]
            for cell in V.nodes[i]["cells"]:
                if cell in V.nodes[n]["cells"]:
                    cells.append(cell)
                else:
                    pass
            
            c1,c2=cp.get_centroid(V, C, cells[0]),cp.get_centroid(V, C, cells[1])
            ### now define a and b
            a=np.array(r)-np.array(c1)
            b=np.array(r)-np.array(c2)
            ma=distance.euclidean(a,[0,0,0])
            mb=distance.euclidean(b,[0,0,0])
            ua=(1/ma)*a
            ub=(1/mb)*b
            
            cos_theta=np.dot(ua,ub)
            
            force=(ua*(cos_theta-1/mb)+ub*(cos_theta-1/ma))*k
            tot_force=tot_force+force 
        
    else:
        if len(list(V[i]))==3:
            r= 0.5*np.array(cp.posi(V,i))+0.5*np.array(cp.posi(V,n))
            cells=[]
            for cell in V.nodes[i]["cells"]:
                if cell in V.nodes[n]["cells"]:
                    cells.append(cell)
                else:
                    pass
            
            c1,c2=cp.get_centroid(V, C, cells[0]),cp.get_centroid(V, C, cells[1])
            ### now define a and b
            a=np.array(r)-np.array(c1)
            b=np.array(r)-np.array(c2)
            ma=distance.euclidean(a,[0,0,0])
            mb=distance.euclidean(b,[0,0,0])
            ua=(1/ma)*a
            ub=(1/mb)*b
            
            cos_theta=np.dot(ua,ub)
            
            force=(ua*(cos_theta-1/mb)+ub*(cos_theta-1/ma))*k
            tot_force=tot_force+force 
    
    return tot_force
    


def bending_force2(V,C,i,k):
    if V.nodes[i]["anchor"]==False:
        r=np.array(cp.posi(V,i))
        
        cell_list=V.nodes[i]["cells"]
        a,b,c=cp.get_centroid(V,C,cell_list[0]),cp.get_centroid(V,C,cell_list[1]),cp.get_centroid(V,C,cell_list[2])
        
        center=(1/3)*(a+b+c)
        
        u=np.array(a)-np.array(b)
        v=np.array(b)-np.array(c)
        mu,mv=np.sqrt(u.dot(u)),np.sqrt(v.dot(v))
        n=np.cross(u,v)/(mu*mv)
        mn=np.sqrt(n.dot(n))
        n=n/mn
        
        if np.dot(n,r-a)>0:
            n=-n
        else:
            pass
        ### unit normal points away from r
        
        l=np.abs(np.dot(r-a,n))*np.array(n)
        
        ## get triangle area
        
        area=(1/2)*distance.euclidean((np.cross(u,v)),[0,0,0])
        if distance.euclidean(l,[0,0,0])>0.1:
            
            force=k*(1/np.sqrt(area))*l
        else:
            force=np.array([0,0,0])
        force=k*(1/np.sqrt(area))*l
    else:
        force=np.array([0,0,0])
    
    return force
    
def bending_force3(V,C,i,k):
    if V.nodes[i]["anchor"]==False:
        r=np.array(cp.posi(V,i))
        v_list=list(V[i])
        
        #a,b,c=cp.get_centroid(V,C,cell_list[0]),cp.get_centroid(V,C,cell_list[1]),cp.get_centroid(V,C,cell_list[2])
        a,b,c=cp.posi(V,v_list[0]),cp.posi(V,v_list[1]),cp.posi(V,v_list[2])
        center=(1/3)*(np.array(a)+np.array(b)+np.array(c))
        
        u=np.array(a)-np.array(b)
        v=np.array(b)-np.array(c)
        mu,mv=np.sqrt(u.dot(u)),np.sqrt(v.dot(v))
        n=np.cross(u,v)/(mu*mv)
        mn=np.sqrt(n.dot(n))
        n=n/mn
        
        if np.dot(n,r-a)>0:
            n=-n
        else:
            pass
        ### unit normal points away from r
        
        l=np.abs(np.dot(r-a,n))*np.array(n)
        
        ## get triangle area
        
        area=(1/2)*distance.euclidean((np.cross(u,v)),[0,0,0])
        if distance.euclidean(l,[0,0,0])>0.1:
            
            force=k*(1/np.sqrt(area))*l
        else:
            force=np.array([0,0,0])
        
    else:
        if len(list(V[i]))==3:
            r=np.array(cp.posi(V,i))
            v_list=list(V[i])
            
            #a,b,c=cp.get_centroid(V,C,cell_list[0]),cp.get_centroid(V,C,cell_list[1]),cp.get_centroid(V,C,cell_list[2])
            a,b,c=cp.posi(V,v_list[0]),cp.posi(V,v_list[1]),cp.posi(V,v_list[2])
            center=(1/3)*(np.array(a)+np.array(b)+np.array(c))
            
            u=np.array(a)-np.array(b)
            v=np.array(b)-np.array(c)
            mu,mv=np.sqrt(u.dot(u)),np.sqrt(v.dot(v))
            n=np.cross(u,v)/(mu*mv)
            mn=np.sqrt(n.dot(n))
            n=n/mn
            
            if np.dot(n,r-a)>0:
                n=-n
            else:
                pass
            ### unit normal points away from r
            
            l=np.abs(np.dot(r-a,n))*np.array(n)
            
            ## get triangle area
            
            area=(1/2)*distance.euclidean((np.cross(u,v)),[0,0,0])
            if distance.euclidean(l,[0,0,0])>0.1:
                
                force=k*(1/np.sqrt(area))*l
            else:
                force=np.array([0,0,0])
        
        else:
            
            force=np.array([0,0,0])
    
    return force
    
def bending_force_test(V,C,i,k):
    if V.nodes[i]["anchor"]==False:
        r=np.array(cp.posi(V,i))
        
        cell_list=V.nodes[i]["cells"]
        a,b,c=cp.get_centroid(V,C,cell_list[0]),cp.get_centroid(V,C,cell_list[1]),cp.get_centroid(V,C,cell_list[2])
        
        center=(1/3)*(a+b+c)
        
        u=np.array(a)-np.array(b)
        v=np.array(b)-np.array(c)
        mu,mv=np.sqrt(u.dot(u)),np.sqrt(v.dot(v))
        n=np.cross(u,v)/(mu*mv)
        mn=np.sqrt(n.dot(n))
        n=n/mn
        if np.dot(n,r-a)>0:
            n=-n
        else:
            pass
        ### unit normal points away from r
        
        l=np.abs(np.dot(r-a,n))*np.array(n)
        
        ## get triangle area
        
        area=(1/2)*distance.euclidean((np.cross(u,v)),[0,0,0])
        
        force=k*(1/np.sqrt(area))*l
    else:
        force=np.array([0,0,0])
    
    return force,n,l,area
        
  
def twist(V,C,i):
    
      if V.nodes[i]["type"]=="avc":
          pass
      else:
          pz=pos[i][2]
          yd=distance.euclidean(mid,pz)
          twist_factor=np.sin(np.pi*((yd-4*a)/mid))
          if twist_factor>0:
              x_array=[]
              y_array=[]
              for p1 in V:
                  if pz-a/2<pos[p1][2]<pz+a/2:
                      x_array.append(pos[p1][0])
                      y_array.append(pos[p1][0])
                  else:
                      pass
                  
              xc,yc=np.mean(x_array),np.mean(y_array)
             
             
              
                 
                 
             
            
             
             
              if  posi(V,i)[2]>mid:
                  
                  
                  ##get angle
                  x,y,z=posi(V,i)[0]-xc,posi(V,i)[1]-yc,posi(V,i)[2]
                  radius=np.sqrt(x**2+y**2)
                  ### set angle correctly
                  if y>0:
                      angle=np.arccos(x/radius)
                  else:
                      angle=2*np.pi-np.arccos(x/radius)
                      
                  new_angle=angle+twist_factor*np.pi/(300*3)
                  
                  x1,y1=radius*np.cos(new_angle)+xc,radius*np.sin(new_angle)+yc
                  V.nodes[i]["pos"]=np.array([x1,y1,z])
              else:
                  x,y,z=posi(V,i)[0]-xc,posi(V,i)[1]-yc,posi(V,i)[2]
                  radius=np.sqrt(x**2+y**2)
                  ### set angle correctly
                  if y>0:
                      angle=np.arccos(x/radius)
                  else:
                      angle=2*np.pi-np.arccos(x/radius)
                      
                  new_angle=angle-twist_factor*np.pi/(300*3)
                  
                  x1,y1=radius*np.cos(new_angle)+xc,radius*np.sin(new_angle)+yc
                  V.nodes[i]["pos"]=np.array([x1,y1,z])
              

def slm_force(V,C,vertex,K):
    rm=posi(V,vertex)
    n_list=list(V[vertex])
    total_force=np.array([0,0,0])
    for i in n_list:
        rn=posi(V,i)
        r_unit=-cp.unit_vector(V, i, vertex)
        r_mag=distance.euclidean(rm,rn)
        force=K*np.array(r_unit)*(r_mag-a/4)
        total_force=total_force+force
    
    return total_force
        

def slm_force_anchor(V,vertex,k):
    rm=posi(V,vertex)
    rn=np.array(V.nodes[vertex]["anchor_pos"])
    
    
    
            
        
        
        
    r_mag=distance.euclidean(rm,rn)
    
    if r_mag>0.01:
        
        r_unit=(1/r_mag)*(np.array(rn)-np.array(rm))
        # if V.nodes[vertex]["type"]=="avc":
        #     force=k*np.array(r_unit)*(r_mag-l/2)
        # else:
        force=k*np.array(r_unit)*(r_mag)
    else:
        force=np.array([0,0,0])
        
        
    
        
    return force