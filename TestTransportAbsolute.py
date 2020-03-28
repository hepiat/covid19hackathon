# -*- coding: utf-8 -*-
"""
Network transportation of the virus - testing with a simple network and the 
corresponding transport through the explicit exchange of people with constant
amount of people

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



# definition of the number of people in each cell
people = np.random.random(len(ind)) * 10000



#%% define the people exchange matrix
# define the amount of people exchanging
exchange_amount = 0.1

# define the matrix
ExchangeMatrix = np.zeros(dist.shape)
for ii in range(len(ind)):
    
    ExchangeMatrix[:,ii] = exchange_amount * np.minimum(np.repeat(people[ii],len(ind)), people)

# limite the transport matrix based on the distance matrix
threshold = 0.3
DiffusionMatrix = 1 / (dist+1)
DiffusionMatrix[DiffusionMatrix < threshold] = 0
DiffusionMatrix_normalized = np.zeros(DiffusionMatrix.shape)
for ii in range(len(ind)):
    
    DiffusionMatrix_normalized[:,ii] =  DiffusionMatrix[:,ii] / np.sum(DiffusionMatrix[:,ii])
    
# limit the exchange matrix by multiplying with the diffusion matrix
ExchangeMatrix_modified = ExchangeMatrix * (DiffusionMatrix_normalized > 0)

# normalize so that the i,i element represents the number of people leaving
for ii in range(len(ind)):
    ExchangeMatrix_modified[ii,ii] = ExchangeMatrix_modified[ii,ii] - np.sum(ExchangeMatrix_modified[:,ii])



#%% def: the update of the state
def StateUpdate(StateOld, ExchangeMatrix_modified):
    
    # define how many people are moving around
    StateRel = np.zeros(StateOld.shape)
    for ii in range(StateRel.shape[0]):
        StateRel[ii,:] = StateOld[ii,:] / np.sum(StateOld[ii,:])
    dState = np.zeros(StateOld.shape)
    for ii in range(StateOld.shape[1]):
        dState[:,ii] = np.matmul(ExchangeMatrix_modified, StateRel[:,ii])
    
    # define the new state
    StateNew = StateOld + dState 
    
    # return the new state
    return StateNew



#%% test the diffusion
# number of time steps
Ttotal = 10

# define the initial state
State = np.zeros((len(ind),2))
State[:,0] = 1
State[0,:] = (0,1)
State = State * np.reshape(np.repeat(people, 2), (len(ind),2))

# store the states over time
StateAll = np.zeros((len(ind), 2, Ttotal))

# loop over the time steps
StateOld = State
for ii in range(Ttotal):
    
    StateNew = StateUpdate(StateOld, ExchangeMatrix_modified)
    
    StateAll[:,:,ii] = StateNew
    
    StateOld = StateNew
    


#%% plot the diffusion of the state 1 people
plt.figure()
for ii in range(Ttotal):
    
    plt.imshow(np.reshape(StateAll[:,1,ii], X.shape), vmax=1000, vmin=0)
    plt.savefig('images_absolute\\' + str(ii) + '.png')
    plt.close()

    

#%% plot also the relative number of state 1 people
plt.figure()
for ii in range(Ttotal):
    
    plt.imshow(np.reshape(StateAll[:,1,ii] / np.sum(StateAll[:,:,ii],axis=1), X.shape), vmax=0.1, vmin=0)
    plt.savefig('images_absolute\\rel_' + str(ii) + '.png')
    plt.close()