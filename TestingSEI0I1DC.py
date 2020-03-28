# -*- coding: utf-8 -*-
"""
computing SEIR results for a single population with the following extensions
- infectious people go through 2 phases: being infectious but not knowing it (I0)
        and being infectious and knowing it (I1)
- removed are split into cured (C) and dead (D)

Created on Fri Mar 27 21:28:54 2020

@author: Roland
"""



#%% libraries
import numpy as np
import matplotlib.pyplot as plt



#%% define the SIR model
def SEI0I1DCmodel(S,E,I0,I1,C,D, beta0, beta1, gamma0, gamma1, eps1, delta, a):
    
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
    
    # return the updates states
    return S, E, I0, I1, C, D
    


#%% create the model properties
# beta: amount of susceptible getting exposed
beta0 = 0.2
beta1 = 0.2

# gamma: amount of infected getting cured
gamma0 = 0.1
gamma1 = 0.1

# delta: people moving from I0 to I1
delta = 0.5

# eps: people dying from I1
eps1 = 0.02

# a: incubation parameter with exponential distribution with parameter a ( the 
# incubation period is a^-1) --> exposed to I0
a = 0.1

# define the initial states of the model
S = 99.0
E = 0.0
I0 = 1.0
I1 = 0.0
C = 0.0
D = 0.0

# number of time steps considered
Ttotal = 200



#%% simulate the process over x time steps
# initialize all states over time
Sall = np.zeros((Ttotal))
Eall = np.zeros((Ttotal))
I0all = np.zeros((Ttotal))
I1all = np.zeros((Ttotal))
Call = np.zeros((Ttotal))
Dall = np.zeros((Ttotal))

# define the starting states
Sold = S
Eold = E
I0old = I0
I1old = I1
Cold = C
Dold = D

# loop over the times
for ii in range(Ttotal):
    
    print('iteration ' + str(ii))
    
    # update the states
    Snew, Enew, I0new, I1new, Cnew, Dnew = SEI0I1DCmodel(Sold,Eold,I0old,I1old,Cold,Dold, beta0, beta1, gamma0, gamma1, eps1, delta, a)
    
    # put into the collection of the results
    Sall[ii] = Snew
    Eall[ii] = Enew
    I0all[ii] = I0new
    I1all[ii] = I1new
    Call[ii] = Cnew
    Dall[ii] = Dnew
    
    # overwrite the old values for the next iteration
    Sold = Snew
    Eold = Enew
    I0old = I0new
    I1old = I1new
    Cold = Cnew
    Dold = Dnew
    
  
    
#%% plot the progress
plt.figure()
plt.plot(Sall, 'g')
plt.plot(Eall, 'y')
plt.plot(I0all, 'r')
plt.plot(I1all, 'm')
plt.plot(Call, 'b')
plt.plot(Dall, 'k')
plt.plot(Sall+Eall+I0all+I1all+Call+Dall, 'g--')
plt.yscale('log')
plt.title('Percentage in log scale')
plt.legend(('susceptible', 'exposed', 'silent infected', 'infected', 'cured', 'dead'))


#%% plot also the change of the state variables
plt.figure()
plt.plot(np.diff(Sall), 'g')
plt.plot(np.diff(Eall), 'y')
plt.plot(np.diff(I0all), 'r')
plt.plot(np.diff(I1all), 'm')
plt.plot(np.diff(Call), 'b')
plt.plot(np.diff(Dall), 'k')
plt.yscale('log')
plt.title('Increase in log scale')
plt.legend(('susceptible', 'exposed', 'silent infected', 'infected', 'cured', 'dead'))