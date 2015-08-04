import numpy as np
import scipy as sp
import pylab as pl
import scipy.stats as stats

import sys

sys.path.append( '../../Exercise_7/source' )

from draw_map import *


def electrostatic_potential(coeffPath):
	#First we need to set up a grid over the high latitudes
	lamb_start = 60
	lamb_end = 90
	phi_start = 0
	phi_end = 361
	N = 8

	grid = np.zeros([lamb_end - lamb_start, phi_end - phi_start])	#Spatial grid
	legendre = np.zeros((30,N + 1,N +1))							#Legendre polynomials (30 arguments, 8 l, and 8 m)
	A = np.zeros((N+1, N+1))										#Coeff arrays A[l,m]
	B = np.zeros((N+1, N+1))

	#Now we need to read in the harmonic expansion coefficients
	#If the coefficient is for A or B is stored in the coeff_type list
	#	l 	m	coeff
	#	#	#	#####
	# The first column gives the l_number (degree), the second the m number (order) and the third is the actual coefficient
	coeff_type = np.genfromtxt(coeffPath, dtype = str, skiprows = 1,  usecols = (0))
	l_number = np.genfromtxt(coeffPath, dtype = int, skiprows = 1,  usecols = (1))
	m_number = np.genfromtxt(coeffPath, dtype = int, skiprows = 1,  usecols = (2))
	coefficients = np.genfromtxt(coeffPath, skiprows = 1,  usecols = (3))

	nCoeff = coefficients.shape[0]

	#Fill the coefficients into the A and B arrays
	for l in range(0,N+1):
			for m in range(0,l+1):
				#Now that we are cycling through the double sum, we need to find the A_lm and B_lm out of the arrays
				#An improvement can be made by using the order the coefficients are in instead of if testing (won't bother since it's fast anyway)
				for n in range(0,nCoeff):
					if coeff_type[n] == 'A' and l == l_number[n] and m == m_number[n]:
						# print 'A' + str(l) + str(m)
						A[l,m] = coefficients[n]
					elif coeff_type[n] == 'B' and l == l_number[n] and m == m_number[n]:
						# print 'B' + str(l) + str(m)
						B[l,m] = coefficients[n]


	#We need a function that calculates the cos(theta), as the arguement taken by the Legendre Polynomial
	def cos_theta(lamb):
		return 2.*(lamb - 75.)/30.

	#Fill the Legendre grid with the arguments
	for i in range(legendre.shape[0]):
		legendre[i,:,:] = cos_theta(i + 60)


	#Then we calculate the actual legendre polynomials	
	for i in range(legendre.shape[0]):
		argument = legendre[i,0,0]
		for l in range(N + 1):
			for m in range(N + 1):
				legendre[i,l,m] = sp.special.lpmv(l,m,argument)	#The scipy function takes as input lpmv( order  , degree , argument ). # This seems to be correct instead lpmv(degree, order, argument)
																# http://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.special.lpmv.html
																# https://github.com/scipy/scipy/blob/master/scipy/special/specfun/specfun.f#L7965
																# These two disagree, assuming that they are representing the same function


	# Then we have all the tools to fill the latitude and longitude grid with the potential field
	#	-	1	.. 	359
	#	60
	#	.
	#	.
	#	89

	def Phi(phi,latitude_index, legendre):
		#This function takes as input an longitude phi, an index corresponding to a latitude
		#and an array of legendre(argument, degree, order) polynomials.
		#Then it calculates a double sum 
		# sum_l sum_m [A(l,m)cos(m*phi) + B(l,m)cos(m*phi)]P(l,m,argument)

		#Convert to radians
		phi = phi*np.pi/180
		total_sum = 0


		for l in range(N + 1):
			for m in range(l + 1):
				total_sum += (A[l,m]*np.cos(m*phi) + B[l,m]*np.sin(m*phi))* legendre[latitude_index, l, m]

		return total_sum


	for latitude_index in range(grid.shape[0]):
		for longitude in range(grid.shape[1]):
			grid[latitude_index, longitude] = Phi(longitude, latitude_index, legendre)

	return grid
		
def plot_potential(grid, name):
	#Variables for plotting purposes
	x = np.arange(0, 361)
	y = np.arange(60, 90)
	X , Y = np.meshgrid(x, y)

	#Plotting a contour plot of the electrostatic potential
	pl.figure()
	contourPlot = pl.contourf(X,Y,grid)
	pl.title('Electrostatic potential in the higher latitudes')

	cbar = pl.colorbar(contourPlot)
	pl.xlabel('Longitude $ [ ^\circ] $')
	pl.ylabel('Latitude  $ [ ^\circ] $')

	#Hardcoded the labels
	cbar_labels = np.arange(-24 , 32 + 1, 24-16)

	cbar.ax.set_yticklabels(cbar_labels)
	cbar.set_label(' $\Phi$  [kV]')#, rotation = 0)
	pl.savefig(name)

	#Drawing on a orthogonal projection of the north pole instead
	fig = plt.figure()

	my_map = draw_map()
	longitude, latitude = my_map(X, Y)
	contour = my_map.contourf(longitude, latitude , grid)
	pl.title('Electrostatic Potential')
	cbar = pl.colorbar(contour, orientation='vertical')
	cbar.set_label(' $\Phi$  [V]')#, rotation = 0)
	pl.savefig('map_' + name)



	pl.show()

if __name__ == "__main__":
	figname = 'potential.eps'
	coeffPath = 'heppner_coeffs.txt'
	potential = electrostatic_potential(coeffPath)
	
	plot_potential( potential ,figname)