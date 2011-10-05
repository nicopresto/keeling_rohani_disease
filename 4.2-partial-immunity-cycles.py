#!/usr/bin/env python

####################################################################
###    This is the PYTHON version of program 4.2 from page 123 of  #
### "Modeling Infectious Disease in humans and animals"            #
### by Keeling & Rohani.										   #
###																   #
### It is a N strain SIR disease model, where the strains are      #
### arranged in a circle and each strain offers partial immunity   #
### (in terms of reduced transmission) to its neighbours.		   #
####################################################################

###################################
### Written by Ilias Soumpasis    #
### ilias.soumpasis@ucd.ie (work) #
### ilias.soumpasis@gmail.com	  #
###################################

import scipy.integrate as spi
import numpy as np
import pylab as pl
from matplotlib.font_manager import FontProperties as fmp

N=4;
beta=40.0;
gamma=9.98;
mu=0.02;
a=0.4;
S0=np.array([0.08, 0.1, 0.1, 0.11]);
P0=np.array([0.4, 0.3, 0.3, 0.29]);
lambda0=np.array([0.15, 0.02, 0.1, 0.01]);
ND=MaxTime=200;
TS=1.0

#####################################################################################
### To be compatible with other versions of programs the 
### following options are available. To try some of them
### uncomment the code (remove '#'):
#####################################################################################
### You may also wish to try:
#(N,beta,gamma,mu,a,S0,P0,lambda0,ND)=(4,40,9.98,0.02,0.25,([0.25, 0.14, 0.25, 0.14]),\
#([0.016, 0.55, 0.016, 0.55]),([0.07, 1e-12, 0.07, 1e-12]),200)
INPUT=np.hstack((S0,P0,lambda0))


def diff_eqs(INP,t):  
	'''The main set of equations'''
	Y=np.zeros((3*N))
	V = INP   
	for i in range(N):
		r=np.mod(i+1,N)
		l=np.mod(i+N-1,N)
		Y[i] = mu - V[i]*(V[(2*N)+i]+V[(2*N)+l]+V[(2*N)+r]) - mu*V[i]
		Y[N+i] =  V[i] * (V[(2*N)+l] + V[(2*N)+r]) - V[N+i] * V[(2*N)+i] - mu*V[N+i]
		Y[(2*N)+i] = beta * (V[i] + a * V[N+i]) * V[(2*N)+i] - gamma*V[(2*N)+i] - mu*V[(2*N)+i]	
	return Y   # For odeint

t_start = 0.0; t_end = ND; t_inc = TS
t_range = np.arange(t_start, t_end+t_inc, t_inc)
RES = spi.odeint(diff_eqs,INPUT,t_range)

print RES

#Ploting
pl.subplot(311)
pl.plot(RES[:,0], 'b-', label='S1')
pl.plot(RES[:,1], 'g-', label='S2')
pl.plot(RES[:,2], 'r-', label='S3')
pl.plot(RES[:,3], 'c-', label='S4')
#pl.xlabel ('Time')
pl.ylabel ('Susceptible')
pl.legend(loc=1, prop=fmp(size='smaller'))
pl.subplot(312)
pl.plot(RES[:,N+0], 'b-', label='P1')
pl.plot(RES[:,N+1], 'g-', label='P2')
pl.plot(RES[:,N+2], 'r-', label='P3')
pl.plot(RES[:,N+3], 'c-', label='P4')
#pl.xlabel ('Time')
pl.legend(loc=1, prop=fmp(size='smaller'))
pl.subplot(313)
pl.plot(RES[:,2*N+0], 'b-', label='l1')
pl.plot(RES[:,2*N+1], 'g-', label='l2')
pl.plot(RES[:,2*N+2], 'r-', label='l3')
pl.plot(RES[:,2*N+3], 'c-', label='l4')
pl.xlabel ('Time')
pl.ylabel ('Force of Infection')
pl.legend(loc=1, prop=fmp(size='smaller'))
pl.show()