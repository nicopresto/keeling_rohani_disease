#!/usr/bin/env python

####################################################################
###    This is the PYTHON version of program 3.1 from page 58 of   #
### "Modeling Infectious Disease in humans and animals"            #
### by Keeling & Rohani.					   #
### It is the SIS model with two different risk-groups.		   #
####################################################################

###################################
### Written by Ilias Soumpasis    #
### ilias.soumpasis@ucd.ie (work) #
### ilias.soumpasis@gmail.com	  #
###################################

import scipy.integrate as spi
import numpy as np
import pylab as pl

beta=[10., 0.1, 0.1, 1.]
gamma=1.0
nH=0.2
IH=1e-5
IL=1e-3
nT=1.0
nL=nT-nH
SH=nH-IH
SL=nL-IL
ND=15.
TS=1.0

INPUT = (SH,IH,SL,IL)

def diff_eqs(INP,t):  
	'''The main set of equations'''
	Y=np.zeros((4))
	V = INP   
	Y[0] = - (beta[0] * V[1] + beta[1] * V[3]) * V[0] + gamma * V[1]
	Y[1] = (beta[0] * V[1] + beta[1] * V[3]) * V[0] - gamma * V[1]
	Y[2] = - (beta[2] * V[1] + beta[3] * V[3]) * V[2] + gamma * V[3]
	Y[3] = (beta[2] * V[1] + beta[3] * V[3]) * V[2] - gamma * V[3]
	return Y   # For odeint

t_start = 0.0; t_end = ND; t_inc = TS
t_range = np.arange(t_start, t_end+t_inc, t_inc)
RES = spi.odeint(diff_eqs,INPUT,t_range)

print RES

#Ploting
pl.subplot(211)
pl.plot(RES[:,1], '-r', label='High Risk')
pl.plot(RES[:,3], '--r', label='Low Risk')
pl.legend(loc=0)
pl.title('Program_3_1.py')
pl.xlabel('Time')
pl.ylabel('Infectious')
pl.subplot(212)
pl.semilogy(RES[:,1], '-r', label='High Risk')
pl.semilogy(RES[:,3], '--r', label='Low Risk')
pl.legend(loc=0)
pl.xlabel('Time')
pl.ylabel('Infectious')
pl.show()
