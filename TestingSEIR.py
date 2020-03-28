# -*- coding: utf-8 -*-
"""
computing SEIR results for a single population 

Created on Fri Mar 27 21:28:54 2020

@author: Roland
"""



#%% libraries
import numpy as np
import matplotlib.pyplot as plt



#%% define the SIR model
def SEIRmodel(S,E,I,R, beta, gamma, a):
    
    N = S + E + I + R
    
    # change of the susecptible
    dSdt = -beta * I * S / N 
    
    # chance of the exposed 
    dEdt = beta * I * S / N - a * E 
    
    # change of the infected
    dIdt = a * E - gamma * I
    
    # change of the removed
    dRdt = gamma * I
    
    # update the states
    S = S + dSdt
    E = E + dEdt
    I = I + dIdt
    R = R + dRdt
    
    # return the updates states
    return S, E, I, R
    


#%% create the model properties
# beta: amount of suscebtible gettting infected 
beta = 0.1

# gamma: amount of infected getting removed
gamma = 0.01

# a: incubation parameter with exponential distribution with parameter a ( the 
# incubation period is a^-1) 
a = 0.1

# define the initial states of the model
S = 95
E = 5
I = 0
R = 0

# number of time steps considered
Ttotal = 100



#%% simulate the process over x time steps
# initialize all states over time
Sall = np.zeros((Ttotal))
Eall = np.zeros((Ttotal))
Iall = np.zeros((Ttotal))
Rall = np.zeros((Ttotal))

# define the starting states
Sold = S
Eold = E
Iold = I
Rold = R

# loop over the times
for ii in range(Ttotal):
    
    print('iteration ' + str(ii))
    
    # update the states
    Snew, Enew, Inew, Rnew = SEIRmodel(Sold, Eold, Iold, Rold, beta, gamma, a)
    
    # put into the collection of the results
    Sall[ii] = Snew
    Eall[ii] = Enew
    Iall[ii] = Inew
    Rall[ii] = Rnew
    
    # overwrite the old values for the next iteration
    Sold = Snew
    Eold = Enew
    Iold = Inew
    Rold = Rnew
    
  
    
#%% plot the progress
plt.figure()
plt.plot(Sall, 'g')
plt.plot(Eall, 'y')
plt.plot(Iall, 'r')
plt.plot(Rall, 'k')