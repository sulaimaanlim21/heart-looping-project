# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 18:15:35 2021

@author: sulai
"""





import numpy
from matplotlib.pyplot import pause
import pylab
import imageio



import matplotlib.animation as animation

import os
import glob
from PIL import Image

from moviepy.editor import *

# files = ['1.png', '2.png', '3.png', '4.png']
# clip = ImageSequenceClip(files, fps = 4) 
# clip.write_videofile("video.mp4", fps = 24)
#image_path = "C:\Users\OneDrive\Documents\PhD\jan\hex\3d\3d_alt\movies"

#os.mkdir(image_path)

path='C:/Users/sulai/OneDrive - Imperial College London/PhD/Simulation/big_files/'
file_name = 'stage_new'

images=[]
step=0
for i in range(200):
    #file_name = 'tn' + str(i+299)
    image=path+file_name + str(i)+'.png'
    images.append(image)
    step+=1
    
clip=ImageSequenceClip(images,fps=30)
clip.write_videofile("C:/Users/sulai/OneDrive - Imperial College London/PhD/Simulation/movies/movie_new.mp4", fps=20)

