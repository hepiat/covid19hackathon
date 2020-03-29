# -*- coding: utf-8 -*-
"""
Network transportation of the virus - testing with a simple network and the 
corresponding transport / diffusion matrix

Created on Sat Mar 28 15:27:37 2020

@author: Roland
"""


#%% import the libraries
import numpy as np
import matplotlib.pyplot as plt



#%% define the geometry 
# the mesh and the corresponding coordinates
x = np.arange(10)
X,Y = np.meshgrid(x,x)

# list the locations
x = np.concatenate(X)
y = np.concatenate(Y)
ind = np.arange(len(x))

# get the distance between the different locations
dist = np.zeros((len(ind), len(ind)))
for ii in range(len(ind)):
        
    dist[ii,:] = np.sqrt((x[ii]-x)**2 + (y[ii]-y)**2)



#%% define the diffusion matrix
# define the threshold for zero diffusion
threshold = 0.3

# construct the matrix
DiffusionMatrix = 1 / (dist+1)
DiffusionMatrix[DiffusionMatrix < threshold] = 0

# normalize the matrix
DiffusionMatrix_normalized = np.zeros(DiffusionMatrix.shape)
for ii in range(len(ind)):
    
    DiffusionMatrix_normalized[:,ii] =  DiffusionMatrix[:,ii] / np.sum(DiffusionMatrix[:,ii])
    


#%% test the diffusion
# number of time steps
Ttotal = 1000

# define the initial state
State = np.zeros(len(ind))
State[0] = 1

# store the states over time
StateAll = np.zeros((len(ind), Ttotal))

# loop over the time steps
StateOld = State
for ii in range(Ttotal):
    
    StateNew = np.matmul(DiffusionMatrix_normalized, StateOld)
    
    StateAll[:,ii] = StateNew
    
    StateOld = StateNew
    


#%% plot the diffusion 
plt.figure()
for ii in range(50): #range(Ttotal):
    
    plt.imshow(np.reshape(StateAll[:,ii], X.shape))
    plt.savefig('images\\' + str(ii) + '.png')
    plt.close()