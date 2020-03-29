# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 16:37:48 2020

@author: Roland
"""

#%% load the libraries 
import numpy as np
import gdal
import matplotlib.pyplot as plt


#%% define the grid
lon = np.arange(-180, 180, 1)
lat = np.arange(-90, 90, 1)
LonGrid, LatGrid = np.meshgrid(lon, lat)

# create a list and indices for all cells
LonList = np.concatenate(LonGrid)
LatList = np.concatenate(LatGrid)
Ncells = len(LatList)
print('Number of cells is ' + str(Ncells))



#%% define the country belonging to the grid points
Country = np.zeros(LonGrid.shape)



#%% define the population map and assigm them to the grid
# load the file
FileName = 'population_density\\gpw_v4_population_density_rev11_2020_1_deg.tif'
ds = gdal.Open(FileName)
geoMatrix = ds.GetGeoTransform()
band = ds.GetRasterBand(1)
arr = band.ReadAsArray()
del ds, band

# get the details of the transform
ulX = geoMatrix[0]
ulY = geoMatrix[3]
xDist = geoMatrix[1]
yDist = geoMatrix[5]
rtnX = geoMatrix[2]
rtnY = geoMatrix[4]
Ny, Nx = arr.shape

# flip the map 
if yDist < 0:
    arr = np.flipud(arr)

# flag for the -inf valued cells for later
FlagInf = arr < -10000

# total population in the world
PopTotal = 8e9

# normalize the population to reflect the world population of 8 billion
Population = arr
Population[FlagInf == 0] = Population[FlagInf == 0] / np.sum(Population[FlagInf == 0]) * PopTotal

# show the population grid for visual check
plt.figure()
plt.imshow(np.log10(Population), vmin = 0, vmax = np.log10(np.max(Population)), origin="lower")



#%% define the health index per grid cell
Health = np.zeros(LatGrid.shape)+1


#%% define the government interventiones per grid cell
Government = np.zeros(LatGrid.shape)


#%% check that the different matrices are consistent
# grid and population matrix
assert np.all(LonGrid.shape == Population.shape), 'inconsistency in between Grid and Population!!'

# ...and the health matrix
assert np.all(LonGrid.shape == Health.shape), 'inconsistency in between Grid and Population!!'

# ...and the government matrix
assert np.all(LonGrid.shape == Government.shape), 'inconsistency in between Grid and Population!!'


#%% remove the inactive cells for efficiency
# based on the information above, remove the cells that do not dbelong to any country and do not have a population
# find cells with no population
FlagNoPopulation = Population <= 0
print('Fraction of cells without any population: ' + str(np.mean(FlagNoPopulation)))

# find the cells with no country attributed
FlagNoCountry = Country < 0
print('Fraction of cells without any attributed country: ' + str(np.mean(FlagNoCountry)))

# combine the 2 flags
FlagRemove = FlagNoPopulation + FlagNoCountry
print('Fraction of cells removed: ' + str(np.mean(FlagRemove >= 1)))
print('Remaining are ' + str(np.sum(FlagRemove >= 1)) + ' cells')

# remove the corresponding cells
LatListClean = LatList[np.concatenate(FlagRemove) >= 1]
LonListClean = LonList[np.concatenate(FlagRemove) >= 1]
PopulationClean = np.concatenate(Population)[np.concatenate(FlagRemove) >= 1]
CountryClean = np.concatenate(Country)[np.concatenate(FlagRemove) >= 1]
NcellsClean = len(CountryClean)



#%% definition of the epidemological model
def SEI0I1DCmodel(S,E,I0,I1,C,D, beta0, beta1, gamma0, gamma1, eps1, delta, a):
    
    # the entire population
    N = S + E + I0 + I1 + C + D
    
    # change of the susecptible
    dSdt = -beta0 * I0 * S / N - beta1 * I1 * S / N 
    
    # chance of the exposed 
    dEdt = beta0 * I0 * S / N  + beta1 * I1 * S / N - a * E 
    
    # change of the infected I0
    dI0dt = a * E - delta * I0 - gamma0 * I0
    
    # change of the infected I1
    dI1dt = delta * I0 - gamma1 * I1 - eps1 * I1
    
    # change of the cured
    dCdt = gamma0 * I0 + gamma1 * I1
    
    # change of the dead
    dDdt = eps1 * I1
    
    # update the states
    S = S + dSdt
    E = E + dEdt
    I0 = I0 + dI0dt
    I1 = I1 + dI1dt
    C = C + dCdt
    D = D + dDdt
    
    # return the updated states
    return S, E, I0, I1, C, D



#%% definition of the transport model
# define the state update as a function of the state StateOld and the ExchangeMatrix
def StateUpdate(StateOld, ExchangeMatrix):
    
    # get the relative distribution of people in each classes 
    StateRel = np.zeros(StateOld.shape)
    for ii in range(StateRel.shape[0]):
        StateRel[ii,:] = StateOld[ii,:] / np.sum(StateOld[ii,:])
    
    # compute the delta state
    dState = np.zeros(StateOld.shape)
    for ii in range(StateOld.shape[1]):
        dState[:,ii] = np.matmul(ExchangeMatrix, StateRel[:,ii])
    
    # define the new state
    StateNew = StateOld + dState 
    
    # return the new state
    return StateNew

# define the ExchangeMatrix describing the number of people transferring from 1 cell to the next in each time step
ExchangeMatrix = np.zeros((NcellsClean, NcellsClean))



#%% simulate the evolution of the virus
# the total number of time steps for the simulation 
Ttotal = 10

# initialize the model parameters 
# beta: amount of susceptible getting exposed
beta0 = 0.2
beta1 = 0.2

# gamme: gamma: amount of infected getting cured
gamma0 = 0.1
gamma1 = 0.1

# eps: people dying from I1
eps1 = 0.01

# delta: people moving from I0 to I1
delta = 0.5

# a: incubation parameter with exponential distribution with parameter a (the incubation period is a^-1) --> exposed to I0
a = 0.1

# and put them into a single dictionary for later
OptsEpid = {'beta0':beta0, 'beta1':beta1, 'gamma0':gamma0, 'gamma1':gamma1, 'gamma0':gamma0, 'delta':delta, 'eps1':eps1, 'a':a}

# and put the transport properties also into a structure
OptsTrans = {'ExchangeMatrix':ExchangeMatrix}

# define the single iteration of the simulation
def SimulateSingleStep(StateOld, OptsEpid, OptsTrans):
    
    # extract the states
    Sold = StateOld[:,0] 
    Eold = StateOld[:,1]
    I0old = StateOld[:,2]
    I1old = StateOld[:,3]
    Cold = StateOld[:,4]
    Dold = StateOld[:,5]
    
    # update the epidemic parameters
    S,E,I0,I1,C,D = SEI0I1DCmodel(Sold,Eold,I0old,I1old,Cold,Dold, 
                                  OptsEpid['beta0'], OptsEpid['beta1'], OptsEpid['gamma0'], OptsEpid['gamma1'], 
                                  OptsEpid['eps1'], OptsEpid['delta'], OptsEpid['a'])
    
    # and put them back into the State structure
    State = np.zeros(StateOld.shape)
    State[:,0] = S
    State[:,1] = E
    State[:,2] = I0
    State[:,3] = I1
    State[:,4] = C
    State[:,5] = D
    
    # update the transport
    StateNew = StateUpdate(State, OptsTrans['ExchangeMatrix'])
    
    return StateNew

# define the initial state for all grid cells
State0 = np.zeros((NcellsClean, 6))
State0[:,0] = PopulationClean

# run the simulation by looping through the simulation steps
# define the storage structure
StateAll = np.zeros((State0.shape[0], State0.shape[1], Ttotal))

# start the iterations
StateOld = State0
for tt in range(Ttotal):
    
    # show the progress
    print('Iterations ' + str(tt) + ' of ' + str(Ttotal) + ' days...')
    
    # update the state
    StateNew = SimulateSingleStep(StateOld, OptsEpid, OptsTrans)
    
    # report on the new state
    StateAll[:,:,tt]
    
    # store the new state
    StateOld = StateNew
    