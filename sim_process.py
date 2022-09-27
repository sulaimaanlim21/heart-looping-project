# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 15:41:49 2022

@author: sulai
"""

import numpy as np
import pickle as pkl

import networkx as nx
from scipy.spatial import distance

import matplotlib.pyplot as plt

from tissue import tissue_generation
import plotly.graph_objects as go
import plotly.express as px


from t1_transition import t1_transition


from sim_stage import sim_stage
from sim_initialise import sim_initialise

def sim_process(shape,parameters,tilt,name):
    V,C = sim_initialise(shape,tilt,name)
    
    V1,C1=sim_stage(V,C,parameters,name)
    
    return V1,C1 