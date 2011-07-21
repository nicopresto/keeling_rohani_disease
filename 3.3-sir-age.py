#!/usr/bin/env python

####################################################################
###    This is the PYTHON version of program 3.3 from page 79 of   #
### "Modeling Infectious Disease in humans and animals"            #
### by Keeling & Rohani.										   #
###																   #
### It is the SIR model with two different age-groups.			   #
####################################################################

###################################
### Written by Ilias Soumpasis    #
### ilias.soumpasis@ucd.ie (work) #
### ilias.soumpasis@gmail.com	  #
###################################

import scipy.integrate as spi
import numpy as np
import pylab as pl

lC=0.066667
mu=[0.0, 0.0166667]
S0=[0.1, 0.1]
I0=[0.0001, 0.0001]
ND=MaxTime=100
beta=np.array([100.0, 10.0, 10.0, 20.0])
gamma=10.0
n0=mu[1]/(lC+mu[1])
n1=1.0-n0
n=np.array([n0,n1])
nu=(lC+mu[1])*n[0]
TS=0.01
INPUT=np.hstack((S0,I0))

def diff_eqs(INP,t):  
	'''The main set of equations'''
	Y=np.zeros((4))
	V = INP    
	Y[0] = nu - (beta[0] * V[2] + beta[2] * V[3]) * V[0] - mu[0] * V[0] - lC * V[0]
	Y[1] = lC * V[0] - (beta[1]*V[2] + beta[3]*V[3]) * V[1] - mu[1] * V[1]
	Y[2] = (beta[0] * V[2] + beta[2] * V[3]) * V[0] - gamma * V[2] - mu[0] * V[2] - lC * V[2]
	Y[3] = (beta[1]*V[2] + beta[3]*V[3]) * V[1] - gamma * V[3] - mu[1] * V[3]
	return Y   # For odeint

t_start = 0.0; t_end = ND; t_inc = TS
t_range = np.arange(t_start, t_end+t_inc, t_inc)
RES = spi.odeint(diff_eqs,INPUT,t_range)

print RES

pl.subplot(211)
pl.plot(RES[:,0], 'g-', label='Ch.Susc')
pl.plot(RES[:,1], 'g--', label='Ad.Susc')
pl.ylabel("Susceptibles")
pl.xlabel('Time')
pl.legend(loc=0)
pl.subplot(212)
pl.plot(RES[:,3], 'r--', label='Ad.Inf')
pl.plot(RES[:,2], 'r-', label='Ch.Inf')
pl.ylabel('Infectious')
pl.xlabel('Time')
pl.legend(loc=0)
pl.show()