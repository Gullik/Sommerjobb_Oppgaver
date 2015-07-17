import numpy as np
import scipy as sp
import pylab as pl
import scipy.stats as stats

from potential_grid import electrostatic_potential

def main():
	#Needed variables
	Re = 6.371E3		#m
	B = np.array([0,0,1]) * 45000.E-9	#tesla

	#Making some needed arrays
	potential = electrostatic_potential()

	#Need 2 grids for the gradient, since it has a 2D direction at each point
	gradient_theta 	= np.zeros(potential.shape)
	gradient_phi	= np.zeros(potential.shape)

	v_theta = np.zeros(potential.shape)
	v_phi = np.zeros(potential.shape)



	#Arrays with the angles, converted to radians for convenience
	theta = np.arange(1,31)
	phi = np.arange(0,360)
	overSin = np.zeros(30)		#To be used in the 1/(rsin(theta)) in the gradient operator for spherical coordinates

	# print theta
	# print phi

	theta = np.deg2rad(theta)
	phi = np.deg2rad(phi)
	sinTheta = np.sin(theta)

	#Need to calculate the differential lengths, need to find out how large a grid step is in meters.
	#For ease we are sticking with 1 grid-node per degree 1 degree
	rho = Re
	dTheta = rho
	dPhi = rho*sinTheta


	#Calculate the gradient for the potential using a centered difference with two gridpoints as the step
	#	The index shifting needed to calculate Phi(i + 1, j) is done for the entire array at the same time 
	# 	[1:-1] 	the inner nodes representing 	Phi(i,j)
	#	[:-2]	nodes 1 node to the left		Phi(i - 1, j)
	#	[2:]	nodes 1 node to the right		Phi(i + 1, j)
	#	And the exact same for the longitude indexes
	# Could be done easier and more foolproof by for-loops, but would be slower

	gradient_theta[1:-1,:] = (potential[2: , :] - potential[:-2 , :])/(rho*2.*dTheta)
	gradient_phi[:,1:-1] = (potential[:,2:] - potential[ :, :-2]    )/(rho*2.*dPhi[:,np.newaxis]*sinTheta[:,np.newaxis])

	# plot_potential_and_electric_field(potential, gradient_theta, gradient_phi)

	#Now we want to calculate the cross products, converting to cartesian coordinates first
	E_x = rho*np.sin(gradient_theta)*np.cos(gradient_phi)
	E_y = rho*np.sin(gradient_theta)*np.sin(gradient_phi)

	v_x = E_y*B[2]
	v_y = -E_x*B[2]

	plot_drift(v_x,v_y)

	pl.show()

def plot_potential_and_electric_field(potential, gradient_theta, gradient_phi):
		#Plotting time...
	#Electrostatic potential
	#Variables for plotting purposes
	x = np.arange(0, 360)
	y = np.arange(60, 90)
	X , Y = np.meshgrid(x, y)
		#Plotting a contour plot of the electrostatic potential
	pl.figure()
	contourPlot = pl.contourf(X,Y,potential)
	pl.title('Electrostatic potential in the higher latitudes')

	cbar = pl.colorbar(contourPlot)
	pl.xlabel('Longitude $ [ ^\circ] $')
	pl.ylabel('Latitude  $ [ ^\circ] $')

	ticks = cbar.ax.get_yticks()
	nTicks = ticks.shape[0]
	cbar_labels = np.arange(-24 , 32 + 1, 24-16)

	cbar.ax.set_yticklabels(cbar_labels)
	cbar.set_label(' $\Phi$  [kV]')#, rotation = 0)
	pl.savefig('Electrostatic_potential')

	#Quiver plot of the electric field E = -grad(Phi), to show the direction in weak areas, all the arrows should have the same length and the magnitude should be given by the color below
	#Getting the magnitude and normalizing
	magnitude = np.sqrt((gradient_theta*gradient_theta) + (gradient_phi*gradient_phi))
	gradient_phi_norm = np.zeros(gradient_phi.shape)
	gradient_phi_norm[1:-1,1:-1] = gradient_phi[1:-1,1:-1]/magnitude[1:-1,1:-1]
	gradient_theta_norm = np.zeros(gradient_theta.shape)
	gradient_theta_norm[1:-1,1:-1] = gradient_theta[1:-1,1:-1]/magnitude[1:-1,1:-1]


	pl.figure()
	x = np.arange(0,potential.shape[1])
	y = 60 + np.arange(0,potential.shape[0])
	X, Y = np.meshgrid(x,y)

	skipArrows = 10
	eArrows = pl.contourf(X,Y, magnitude)
	pl.quiver(X[:,::skipArrows],Y[:,::skipArrows], -gradient_phi_norm[:,::skipArrows], -gradient_theta_norm[:,::skipArrows], angles = 'xy', scale = 40)

	cbar = pl.colorbar(eArrows)

	#Labels
	pl.title('Electric field in the higher latitudes')
	cbar.set_label('$|\\vec{E}|$')
	pl.xlabel('Longitude [ $ ^\circ$]')
	pl.ylabel('Latitude  [ $ ^\circ$]')

	pl.savefig('eArrows.eps')


if __name__ == '__main__':
    main()

