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
    # if C.nodes[cell_index]["type"]=="avc":
    #     gm=1
    # else:
    gm=growth_factor(V,C,cell_index)
    if C.nodes[cell_index]["type"]=="atr":
        A,P=A3,3.4*np.sqrt(A3)
        
    elif C.nodes[cell_index]["type"]=="avc":
        A,P= A1,3.4*np.sqrt(A1)
    else:
        A,P= A2,3.4*np.sqrt(A2)
    area=C.nodes[cell_index]["area"]
    perimeter=C.nodes[cell_index]["perimeter"]
    dA_by_dH=area_deriv(V,C,vertex_index,cell_index)
    dP_by_dH=perimeter_deriv(V,C,vertex_index,cell_index)
    
    dE_by_dH=2*(area-A*gm)*np.array(dA_by_dH) + 2*(perimeter-P)*np.array(dP_by_dH)
    
    return dE_by_dH


def area_deriv_k



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
    if cp.get_centroid(V,C,i)[2]<mid:
        ref_angle=np.pi*(3/2)
    else:
        ref_angle=np.pi*(1/2)
    angle=C.nodes[i]["angle"]   
    if np.pi>angle>0:
        
        growth_factor=1+1/(np.cosh(ref_angle-angle))-growth_z
    else:
        growth_factor=1+1/(np.cosh(ref_angle-angle))+growth_z
    
    # if growth_factor>0:
    growth_modifier=1+growth_factor
    # else:
    #     growth_modifier=1+0.2*growth_factor
    
    return growth_modifier

def growth_factor2(V,C,i):
    if C.nodes[i]["type"]=="atr":
        ref_angle=np.pi*(3/4)
    else:
        ref_angle=np.pi*(7/4)
    angle=C.nodes[i]["angle"]    
    growth_factor=(np.cos(ref_angle-angle)+1)/2
    
    
    
    growth_modifier=1+growth_factor*2/3
    
    
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



def area_perimeter_calculate(V,C):
    for c in C:
        C.nodes[c]["area"]=cp.get_area(V,C,c)
        C.nodes[c]["perimeter"]=cp.get_perimeter(V,C,c)
        


def bending_force_AB(V,C,vertex_index,a,b):
    ### term given by paper 
    
    ### term 1 is (-1/AaAb)[(1/Aa)(dAa/dx)+(1/Ab)(dAb/dx)](Aa.Ab)
    ##term 2 is(1/AaAb)[Aa*dAb*dot+Ab*dAb*dot]
    Aa,Ab=C.nodes[a]["area"],C.nodes[b]["area"]
    dAa,dAb=area_deriv(V, C, vertex_index, a),area_deriv(V, C, vertex_index, b)
    an,bn=C.nodes[a]["normal"],C.nodes[b]["normal"]
    a_dot_b=np.dot(an,bn)
    
    term_1=(-1/(Aa*Ab))*((1/Aa)*dAa+(1/Ab)*dAb)*a_dot_b
    term_2=(1/(Aa*Ab))*(a_dot_b)*(Aa*dAb+Ab*dAa)
    force=term_1+term_2
    return force


       
def bending_force_new(V,C,vertex_index):
    ### term given by paper 
    
    ### term 1 is (-1/AaAb)[(1/Aa)(dAa/dx)+(1/Ab)(dAb/dx)](Aa.Ab)
    
    cell_list=V.nodes[vertex_index]["cells"]
    force=np.array([0,0,0])
    if len(cell_list)==3:
        a,b,c=cell_list[0],cell_list[1],cell_list[2]
        force+=bending_force_AB(V,C,a,b)+bending_force_AB(V,C,c,b)+bending_force_AB(V,C,a,c)
    else:
        a,b=cell_list[0],cell_list[1]
        force+=bending_force_AB(V,C,a,b)
        
    return force
        