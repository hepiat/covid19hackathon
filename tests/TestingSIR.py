# -*- coding: utf-8 -*-
"""
computing SIR results for a single population 

Created on Fri Mar 27 21:28:54 2020

@author: Roland
"""



#%% libraries
import numpy as np
import matplotlib.pyplot as plt



#%% define the SIR model
def SIRmodel(S,I,R, beta, gamma):
    
    N = S + I + R
    
    # change of the susecptible
    dSdt = -beta * I * S / N 
    
    # change of the infected
    dIdt = beta * I * S / N - gamma * I
    
    # change of the removed
    dRdt = gamma * I
    
    # update the states
    S = S + dSdt
    I = I + dIdt
    R = R + dRdt
    
    # return the updates states
    return S, I, R
    


#%% create the model properties
# beta: amount of suscebtible gettting infected 
beta = 0.5

# gamma: amount of infected getting removed
gamma = 0.01

# define the initial states of the model
S = 95
I = 5
R = 0

# number of time steps considered
Ttotal = 100



#%% simulate the process over x time steps
# initialize all states over time
Sall = np.zeros((Ttotal))
Iall = np.zeros(Ttotal)
Rall = np.zeros(Ttotal)

# define the starting states
Sold = S
Iold = I
Rold = R

# loop over the times
for ii in range(Ttotal):
    
    print('iteration ' + str(ii))
    
    # update the states
    Snew, Inew, Rnew = SIRmodel(Sold, Iold, Rold, beta, gamma)
    
    # put into the collection of the results
    Sall[ii] = Snew
    Iall[ii] = Inew
    Rall[ii] = Rnew
    
    # overwrite the old values for the next iteration
    Sold = Snew
    Iold = Inew
    Rold = Rnew
    
  
    
#%% plot the progress
plt.figure()
plt.plot(Sall, 'g')
plt.plot(Iall, 'r')
plt.plot(Rall, 'k')