import numpy as np
import scipy as sp
import pylab as pl
import scipy.stats as stats

from potential_grid import electrostatic_potential

def main():
	#Needed variables
	Re = 6.371E6		#m
	B  = 45000.E-9		#tesla  #Directed in z-direction

	#Making some needed arrays
	potential = electrostatic_potential()

	#Need 2 grids for the gradient, since it has a 2D direction at each point
	gradient_theta 	= np.zeros(potential.shape)
	gradient_phi	= np.zeros(potential.shape)

	#Arrays with the angles, converted to radians for convenience
	theta = np.arange(1,31)
	phi = np.arange(0,360)
	# overSin = np.zeros(30)		#To be used in the 1/(r sin(theta)) in the gradient operator for spherical coordinates

	theta = np.deg2rad(theta)
	phi = np.deg2rad(phi)
	sinTheta = np.sin(theta)

	#Need to calculate the differential lengths, need to find out how large a grid step is in meters.
	#For ease we are sticking with 1 grid-node per degree 1 degree
	rho = Re
	dTheta = 1.
	dPhi = 1.

	#Calculate the gradient for the potential using a centered difference with two gridpoints as the step
	#	The index shifting needed to calculate Phi(i + 1, j) is done for the entire array at the same time 
	# 	[1:-1] 	the inner nodes representing 	Phi(i,j)
	#	[:-2]	nodes 1 node to the left		Phi(i - 1, j)
	#	[2:]	nodes 1 node to the right		Phi(i + 1, j)
	#	And the exact same for the longitude indexes
	# Could be done easier and more foolproof by for-loops, but would be slower

	gradient_theta[1:-1,:] = (potential[2: , :] - potential[:-2 , :])/(rho*2.*dTheta)
	gradient_phi[:,1:-1] = (potential[:,2:] - potential[ :, :-2]    )/(rho * sinTheta[:,np.newaxis] * 2. *dPhi)

	plot_potential_and_electric_field(potential, gradient_theta, gradient_phi)

	print gradient_theta[2,0]

	#Since the electric field is -grad(potential), we use the negative gradient below
	v_theta = cross_with_B(-gradient_theta, theta, phi, B)
	v_phi = cross_with_B(-gradient_phi, theta, phi, B)

	#Stopping here for now
	plot_drift(v_theta,v_phi)


	pl.show()

def cross_with_B(matrix, theta, phi, B):
	#Now we want to calculate the cross products, converting to cartesian coordinates first
	sinTheta = np.sin(theta)
	cosTheta = np.cos(theta)
	sinPhi = np.sin(phi)
	cosPhi = np.cos(phi)

	E_x = matrix * sinTheta[:,np.newaxis] * cosPhi[np.newaxis,:]
	E_y = matrix * sinTheta[:,np.newaxis] * sinPhi[np.newaxis,:]
	E_z = matrix * cosTheta[:,np.newaxis]

	#Do the cross product
	v_x = E_y*B
	v_y = -E_x*B
	v_z = 0


	#Then we return the new velocity magnitude
	return np.sqrt(v_x*v_x + v_y*v_y)/(B*B)




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
	pl.quiver(X[:,::skipArrows],Y[:,::skipArrows], -gradient_phi_norm[:,::skipArrows], -gradient_theta_norm[:,::skipArrows], angles = 'uv', scale = 30)

	cbar = pl.colorbar(eArrows)

	#Labels
	pl.title('Electric field in the higher latitudes')
	cbar.set_label('$|\\vec{E}|$')
	pl.xlabel('Longitude [ $ ^\circ$]')
	pl.ylabel('Latitude  [ $ ^\circ$]')

	pl.savefig('eArrows.eps')

	return

def plot_drift(v_theta, v_phi):

	pl.figure()
	x = np.arange(0,v_theta.shape[1])
	y = 60 + np.arange(0,v_theta.shape[0])
	X, Y = np.meshgrid(x,y)

	skipArrows = 10
	
	Q = pl.quiver(X[:,::skipArrows],Y[:,::skipArrows], v_theta[:,::skipArrows], v_phi[:,::skipArrows], angles = 'uv', scale = 100)
	pl.quiverkey(Q, 1,1.05, 10, '$10$ m/s')

	#Labels
	pl.title('Drift the higher latitudes')
	pl.xlabel('Longitude [ $ ^\circ$]')
	pl.ylabel('Latitude  [ $ ^\circ$]')

	pl.savefig('vArrows.eps')

	return


if __name__ == '__main__':
    main()

