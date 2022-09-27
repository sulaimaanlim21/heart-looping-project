# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 11:56:50 2022

@author: sulai
"""

import numpy as np
import pickle as pkl

import networkx as nx
#from tissue_new import tissue_generation
from sim_initialise import sim_initialise
#from plotting import plotting
from sim_stage_new import sim_stage
from sim_process import sim_process
import cell_properties as cp
from scipy.spatial import distance
from sim_stage_bend import sim_bend
path='C:/Users/sulai/OneDrive - Imperial College London/PhD/Simulation/big_files/'

V,C=nx.read_gpickle(path+'stage_new199.pickle'),nx.read_gpickle(path+'stagec_new199.pickle')

#V,C=sim_initialise([16,16],True,"new")


#V2,C2=sim_bend(V,C,"bend_testy")
V2,C2=sim_stage(V,C,[1.4,1.6,1.6],"new_combined")