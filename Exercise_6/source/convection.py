import numpy as np
import scipy as sp
import pylab as plt
import scipy.stats as stats
import os
import sys

sys.path.append('../../Exercise_5/source')
sys.path.append('../../Exercise_7/source')

from potential_grid import *
from draw_map import draw_map

def convection(potential):
	#Needed variables
	Re = 6.371E6		#m
	B  = 45000.E-9		#tesla  #radially directed



	#Need 2 grids for the gradient, since it has a 2D direction at each point
	gradient_theta 	= np.zeros(potential.shape)
	gradient_phi	= np.zeros(potential.shape)

	#Arrays with the angles, converted to radians for convenience
	theta = np.arange(1,31)
	phi = np.arange(0,360)

	theta = np.deg2rad(31 - theta)		#[60,...90] -> [30,...,0] #Used spherical coordinates
	phi = np.deg2rad(phi)
	sinTheta = np.sin(theta)

	#Need to calculate the differential lengths, need to find out how large a grid step is in meters.
	#For ease we are sticking with 1 grid-node per degree 1 degree
	rho = Re
	dTheta = 1.
	dPhi = 1.

	#Calculate the gradient for the potential using a centered difference with two gridpoints as the step
	#	The index shifting needed to calculate Phi(i + 1, j) is done for the entire array at the same time 
	# 	[1:-1 ,:] 	the inner nodes representing 	Phi(i,j)
	#	[ :-2 ,:]	nodes 1 node to the left		Phi(i - 1, j)
	#	[2:	  ,:]	nodes 1 node to the right		Phi(i + 1, j)
	#	And the exact same for the longitude indexes
	E_theta = np.zeros(potential.shape)
	E_phi = np.zeros(potential.shape)

	E_theta[1:-1,:] = - ( potential[2: , :] - potential[:-2 , :] )/(rho*2.*dTheta)
	E_phi  [:, 1:-1]  = - (potential[:,2:] - potential[ :, :-2]    )/(rho * sinTheta[:,np.newaxis] * 2. *dPhi)

	# plt.t_potential_and_electric_field(potential, E_theta, E_phi)

	#Since the electric field is -grad(potential), we use the negative gradient below
	v_theta, v_phi = cross_with_B(E_theta, E_phi, theta, phi, B)
	v_theta /=(B*B)
	v_phi /= (B*B)

	#Stopping here for now
	# plt.t_drift(v_theta,v_phi)


	return v_theta, v_phi, E_theta, E_phi

def cross_with_B(E_theta, E_phi, theta, phi, B):
	# This takes an latitude and longitude array, along with the latitude and longitude angles as input. 
	# In addition to the magnitude of a radially directed B-field
	# Then it return the cross product as one array of magnitude directed in theta's direction and one in phi's.

	#It can be shown the cross(theta, rho) = - phi evaluated at a point
	#Likewise cross(phi,rho) = theta

	v_theta = -E_phi*B
	v_phi 	= E_theta*B 

	return v_theta, v_phi




def plot_electric_field(E_theta, E_phi):
	#Variables for plotting purposes
	x = np.arange(0, 360)
	y = np.arange(60, 90)
	X , Y = np.meshgrid(x, y)

	#Quiver plot of the electric field E = -grad(Phi), to show the direction in weak areas, all the arrows should have the same length and the magnitude should be given by the color below
	#Getting the magnitude and normalizing
	magnitude = np.sqrt((E_theta*E_theta) + (E_phi*E_phi))
	E_phi_norm = np.zeros(E_phi.shape)
	E_phi_norm[1:-1,1:-1] = E_phi[1:-1,1:-1]/magnitude[1:-1,1:-1]
	E_theta_norm = np.zeros(E_theta.shape)
	E_theta_norm[1:-1,1:-1] = E_theta[1:-1,1:-1]/magnitude[1:-1,1:-1]

	plt.figure()
	x = np.arange(0,potential.shape[1])
	y = 60 + np.arange(0,potential.shape[0])
	X, Y = np.meshgrid(x,y)

	skipArrows = 10
	eArrows = plt.contourf(X,Y, magnitude)
	plt.quiver(X[:,::skipArrows],Y[:,::skipArrows], E_phi_norm[:,::skipArrows], E_theta_norm[:,::skipArrows], angles = 'uv', scale = 40)

	cbar = plt.colorbar(eArrows)

	#Labels
	plt.title('Electric field in the higher latitudes')
	cbar.set_label('$|\\vec{E}|$')
	plt.xlabel('Longitude [ $ ^\circ$]')
	plt.ylabel('Latitude  [ $ ^\circ$]')

	plt.savefig('eArrows.eps')

	#Let's plot the electric field on the map as well
	fig = plt.figure()
	my_map = draw_map()

	longitude, latitude = my_map(X, Y)

	my_map.contourf(longitude, latitude, magnitude)
	my_map.quiver(longitude[:,::skipArrows],latitude[:,::skipArrows], E_phi_norm[:,::skipArrows], E_theta_norm[:,::skipArrows], angles = 'uv', scale = 40)
	plt.savefig('map_efield')

	return 

def plot_drift(v_theta, v_phi):

	# plt.figure()
	x = np.arange(0,v_theta.shape[1])
	y = 60 + np.arange(0,v_theta.shape[0])
	X, Y = np.meshgrid(x,y)

	skipArrows = 10
	
	#Let's normalize the plt.t as well
	plt.figure()

	magnitude = np.sqrt((v_theta*v_theta) + (v_phi*v_phi))
	v_phi_norm = v_phi[1:-1,1:-1]/magnitude[1:-1,1:-1]
	v_theta_norm = v_theta[1:-1,1:-1]/magnitude[1:-1,1:-1]

	eArrows = plt.contourf(X,Y, magnitude)
	Q = plt.quiver(X[:,::skipArrows],Y[:,::skipArrows], v_phi_norm[:,::skipArrows], v_theta_norm[:,::skipArrows], angles= 'uv',  scale = 40)
	cbar = plt.colorbar(eArrows)
	cbar.set_label('$m/s$')

	#Labels
	plt.title('Drift the higher latitudes')
	plt.xlabel('Longitude [ $ ^\circ$]')
	plt.ylabel('Latitude  [ $ ^\circ$]')

	plt.savefig('vArrows.eps')

	#Let us also draw the currents on a map
	fig = plt.figure()
	my_map = draw_map()

	longitude, latitude = my_map(X, Y)

	contour = my_map.contourf(longitude, latitude, magnitude)
	my_map.quiver(longitude[:,::skipArrows],latitude[:,::skipArrows], v_phi_norm[:,::skipArrows], v_theta_norm[:,::skipArrows], angles = 'uv', scale = 40)
	cbar = plt.colorbar(contour)
	cbar.set_label('m/s')

	plt.savefig('map_drift.eps')



	return


if __name__ == '__main__':

	coeffPath = '../../Exercise_5/source/heppner_coeffs.txt'

	#Making some needed arrays
	potential = electrostatic_potential(coeffPath)

	v_theta, v_phi, E_theta, E_phi  = convection(potential)
	# plot_potential(potential, 'potential.eps')
	plot_electric_field(E_theta, E_phi)
	plot_drift(v_theta,v_phi)
	plt.show()

