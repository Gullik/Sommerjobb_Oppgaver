import numpy as np
import scipy as sp
import pylab as pl
import scipy.stats as stats
from time import time

#Define variables
Re = 6.371E3			#m	Earths radius
mu_0 = 4.*np.pi*1.E-7	#N/A^2
resolution = int(1E6)	#Spatial resolution
step = (5*Re - -5*Re)/resolution
omega = 2.*1.E-13		#Hz	Frequency
magConstant = 400.E-9	#T
densityConstant = 50. 	#m^-3

x = np.arange(-5*Re,5*Re, step)

nu = np.zeros(resolution)      #Displacement velocity in y direction
xi = np.zeros(resolution)	   #Displacement
b  = np.zeros(resolution)	   #Magnetic field
rho= np.zeros(resolution)	   #Mass density
va= np.zeros(resolution)	   #Alfven velocity


#Precalculating the magnetic field and density on the spatial nodes used
rho = densityConstant*x**-4
b   = magConstant*x**-3
va 	= b**2/(mu_0*np.abs( rho ))
omega2 = omega*omega
for i in range(0,resolution - 1):	#Need to trim it to avoid problems at 1/x**3 = inf
	if va[i] > 1E-15:
		va[i] = 1E-15


for n in range(0,1):

	print 'Running trial with omega ' + str(omega)
	#Initial conditions
	NU = 100	#displacement velocity 
	XI = 0		#displacement
	DNU = 0		
	DXI = 0     

	#Trial for solving the boundary problem
	for i in range(0,resolution - 1):
		#Change since last step
		DNU = step * (-XI * omega2 / va[i])
		DXI = step * NU

		# print DNU
		#New position and velocity
		XI += DXI
		NU += DNU

		# print NU
		#Storing positions
		xi[i + 1] = XI
		nu[i + 1] = NU


		if i % int(0.1*resolution) == 0:
			print str( float(i) / float(resolution) * 100)  + '%'

	f, ax = pl.subplots(2)
	ax[0].plot(x,xi)
	ax[1].plot(x,nu)

	pl.figure()
	pl.plot(xi,x)

	omega += 0.01

pl.show()
