# This code simulates an ion in bouncing around in the dipolefield aroud the earth 
# using an RK4 method
#
#	python filename.py nSteps timestep
#	
#	E.G
# 	python dipfield_rk4.py 4 -6
#	will start a run with 10^4 steps and 10^-6 timesteps

import sys
import numpy as np
import pylab as plt
from time import time
from mpl_toolkits.mplot3d import Axes3D

def r_length(x,y,z):
	return (x*x + y*y +z*z)**.5

def dipMomz(z):
	return 7.94E22

def bx(x,y,z, C2):
	mz = dipMomz(z)
	return C2*(3.*x*(mz*z))/r_length(x,y,z)**5
	# return 0

def by(x,y,z,C2):
	mz = dipMomz(z)
	return C2*(3.*y*(mz*z))/r_length(x,y,z)**5
	# return 0

def bz(x,y,z,C2):
	mz = dipMomz(z)
	length = r_length(x,y,z)
	return C2*((3.*z*(mz*z))*length**-5 - mz*length**-3)

def b(r,C2):
	#Takes all the magnetic field and returns them as a 3d vector
	x = r[0]
	y = r[1]
	z = r[2]
	 
	return np.array([bx(x,y,z,C2), by(x,y,z,C2), bz(x,y,z,C2)])


def function(y, mag_field, C1, C2):
	#	This function takes as the input a phase vector of positions, and velocities, and a magnetic field.
	# 	Then it calculates the the function f = (v, q/m*vxB) and returns that phase vector
	temp = np.zeros(6)

	temp[:3] = y[3:]														#Calculation of positions dr/dt = v
	temp[3:] = C1*np.cross( y[3:], mag_field(y[3:], C2) )					#Calculation of force 	  dv/dt	= q/m*(v x B)

	return temp


def rk4(y, dydt, C1, C2, steps, timestep):
	#This takes a function and it's derivative and
	#uses the RK4 method, and returns a vector of positions

	k1 = np.zeros(6)
	k2 = np.zeros(6)
	k3 = np.zeros(6)
	k4 = np.zeros(6)

	# for i in range(steps):
	# 	k1 = timestep*dydt( y       , b, C1, C2)
	# 	k2 = timestep*dydt( y + k1/2, b, C1, C2)
	# 	k3 = timestep*dydt( y + k2/2, b, C1, C2)
	# 	k4 = timestep*dydt( y + k3  , b, C1, C2)

	# 	y = timestep/6*( k1 + k2/2 + k3/2 + k4 )

	dydt( y       , b, C1, C2)

	return



if __name__ == '__main__':

	#Defining the physical parameters
	e = 1.602E-19			#C
	Re = 6.371E6			#m
	mu0 = 4.*np.pi*1.E-7	#N/A^2
	q = e
	m = 16*1.674E-26  		#kg

	steps = 10**int(sys.argv[1])
	timestep = 10**int(sys.argv[2])

	C1 = q/m            #Used in the integration
	C2 = mu0/(4*np.pi)	#Constant in the magnetic field calculation
	
	#Setting initial conditions
	r = np.array([5.*Re, 0., 0.]) 		#m   

	pitch 	= np.deg2rad(25)
	v_0 	= 10**7

	v = np.array([v_0*np.cos(pitch), 0.,  v_0*np.sin(pitch)])		#m/s


	y = np.zeros(6)
	y = np.array([r[0], r[1], r[2], v[0], v[1], v[2]])

	dydt = function

	rk4(y, function, C1, C2, steps, timestep)

