import numpy as np
import scipy as sp
import pylab as plt
import scipy.stats as stats
import os
import sys

sys.path.append( '../../Exercise_6/source' )

from convection import *

def field_aligned_currents(v_theta, v_phi):
	#This function takes a input matrices with the latitudal and longitudal velocities,
	#then it returns the field aligned currents, with a magnetic field vertically directed
	#		
	#			FACs = Sigma_P dot(B_0, curl(v))
	#

	Re = 6.371E6		#m
	B0  = 45000.E-9		#tesla  #radially directed
	sigma_P = 1

	theta = np.arange(1,31)
	phi = np.arange(0,360)

	# print theta[2:] - theta[:-2]

	theta = np.deg2rad(31 - theta)		#[60,...90] -> [30,...,0] #Used spherical coordinates
	phi = np.deg2rad(phi)
	sinTheta = np.sin(theta)

	r = Re
	dPhi = 1
	dTheta = 1

	ddTheta = np.zeros(v_theta.shape)
	ddPhi = np.zeros(v_theta.shape)
	curl = np.zeros(v_theta.shape)


	#Computes the derivates with respect to the longitude and latitude by a centered discretization
	# 					d/dx f(x) = [f(x+h) - f(x-h)]/2h
	ddTheta[1:-1,:] = v_phi[2:,:]*sinTheta[2:, np.newaxis] - v_phi[:-2,:]*sinTheta[:-2,np.newaxis]
	ddTheta /= 2.*dTheta
	ddPhi[:,1:-1] 	= v_theta[:,2:] - v_theta[:,:-2]
	ddPhi 	/= 2.*dPhi

	curl = (ddTheta - ddPhi)/(r*sinTheta[:,np.newaxis])

	FACs = - sigma_P * B0 * curl

	return FACs

def plot_FACs(FACs):
	x = np.arange(0, 360)
	y = np.arange(60, 90)

	#Since we the derivatives at the edges hasn't been done properly
	#we'll just shave them of the representation
	FACs[:, 0:2] = 0
	FACs[:, -2:] = 0
	FACs[0:2, :] = 0
	FACs[-2:,  :] = 0

	plt.figure()
	contourPlot = pl.contourf(x,y,FACs)
	plt.title('Field Aligned Currents')

	cbar = pl.colorbar(contourPlot)
	plt.xlabel('Longitude $ [ ^\circ] $')
	plt.ylabel('Latitude  $ [ ^\circ] $')

	cbar.set_label('$ms^{-1}$')

	plt.savefig('facs.eps')

	return

if __name__ == '__main__':
	
	v_theta, v_phi = convection()
	# plot_drift(v_theta,v_phi)
	FACs = field_aligned_currents(v_theta, v_phi)

	plot_FACs(FACs)



	plt.show()

