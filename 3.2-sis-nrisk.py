#!/usr/bin/env python

####################################################################
###    This is the PYTHON version of program 3.2 from page 64 of   #
### "Modeling Infectious Disease in humans and animals"            #
### by Keeling & Rohani.										   #
###																   #
### It is the SIS model with m different risk-groups			   #
### Note we no-longer explicitly model the susceptible class.	   #
####################################################################

###################################
### Written by Ilias Soumpasis    #
### ilias.soumpasis@ucd.ie (work) #
### ilias.soumpasis@gmail.com	  #
###################################

import scipy.integrate as spi
import numpy as np
import pylab as pl

waifw=np.array([0.0, 3.0, 10.0, 60.0, 100.0])
beta=0.0016*np.outer(waifw,waifw)
gamma=0.2*np.array([1.0, 1.0, 1.0, 1.0, 1.0])
n=np.array([0.06, 0.31, 0.52, 0.08, 0.03])
I=np.array([0.0, 0.0, 0.0, 0.0, 1e-5])
ND=30
TS=1.0
m=5
INPUT=I

def diff_eqs(INP,t):  
	'''The main set of equations'''
	Y=np.zeros((m))
	V = INP    
	for i in range(m):
		Y[i] = np.multiply(np.dot(beta[i], V), (n[i]-V[i]))-gamma[i]*V[i]
	return Y   # For odeint

t_start = 0.0; t_end = ND; t_inc = TS
t_range = np.arange(t_start, t_end+t_inc, t_inc)
RES = spi.odeint(diff_eqs,INPUT,t_range)

print RES
pl.subplot(221)
for i in range(m):
	pl.plot(RES[:,i], label = 'risk group %s' %(i+1))
pl.ylabel('Infectious, I')
pl.xlabel('Time')
pl.legend(loc=0)
pl.subplot(222)
for i in range(m-1):
	pl.semilogy(RES[:,i+1], label = 'risk group %s' %(i+2))
pl.ylabel('Infectious, I')
pl.xlabel('Time')
pl.legend(loc=0)
pl.subplot(223)
for i in range(m):
	pl.plot(RES[:,i]/n[i], label = 'risk group %s' %(i+1))
pl.ylabel('Relative level of Infection, I/n')
pl.xlabel('Time')
pl.legend(loc=0)
pl.subplot(224)
for i in range(m-1):
	pl.semilogy(RES[:,i+1]/n[i], label = 'risk group %s' %(i+2))
pl.ylabel('Relative level of Infection, I/n')
pl.xlabel('Time')
pl.legend(loc=0)

pl.show()