import numpy as np
import scipy as sp
import pylab as pl
import scipy.stats as stats
from time import time

##Still needs some work, go through all the constants


#Define variables
Re = 6.371E6			#m	Earths radius
mu_0 = 4.*np.pi*1.E-7	#N/A^2
resolution = int(1E6)	#Spatial resolution
step = (5*Re - (-5*Re))/resolution
omega = 1.E-41			#Hz	Frequency
magConstant = 400.E-9	#T
densityConstant = 50. 	#m^-3





x = np.arange(-5*Re,5*Re, step)

nu = np.zeros(resolution)      #Displacement velocity in y direction
xi = np.zeros(resolution)	   #Displacement
b  = np.zeros(resolution)	   #Magnetic field
rho= np.zeros(resolution)	   #Mass density
vaSquared= np.zeros(resolution)	   #Alfven velocity

trimMag = 1.E-26
trimVA2 = 1.E-60

#Precalculating the magnetic field and density on the spatial nodes used
rho = densityConstant*x*x
b   = magConstant*x**-3				#Could this also be in per cm?
# for i in range(0,resolution):	#Need to trim it to avoid problems at b = 400x**-3 => inf
# 	if np.abs(b[i]) > trimMag:
# 		b[i] = np.sign(b[i]) * trimMag


vaSquared 	= b*b/(mu_0*np.abs( rho ))
# for i in range(0,resolution):	#Need to trim it to avoid problems at 1/x**3 = inf
# 	if vaSquared[i] > trimVA2:
# 		vaSquared[i] = trimVA2


for n in range(0,100):

	print 'Running trial with omega ' + str(omega)
	#Initial conditions
	NU = 100	#displacement velocity 
	XI = 0		#displacement
	DNU = 0		
	DXI = 0

	xi[0] = XI
	nu[0] = NU

	C1 = step * omega*omega / vaSquared   #Precalculating
	# for i in range(0,resolution):	#Need to trim it to avoid problems at 1/x**3 = inf
	# 	if C1[i] > trimNumber:
	# 		C1[i] = trimNumber







	#Trial for solving the boundary problem
	for i in range(0,resolution - 1):
		#Change since last step
		DNU = -XI *C1[i]
		DXI = step * NU

		# if i > 10:
		# 	break

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


	# if np.abs(xi[i]- 0) < 100000:#  and  np.abs(nu[i]) - 100 < 5: 
	f, ax = pl.subplots(2, sharex = True)
	pl.suptitle('Displacement and Spatial Displacement Velocity of a cold plasma')


	ax[0].plot(x,xi)
	ax[1].plot(x,nu)

	# #Change axes
	# xticks = ax[0].get_xticks()
	# locs = ax[0].get_xticks()
	# ax[0].xaxis.set_ticks(locs, map(lambda x: "%.1f" % x, locs*1.E3))

	ax[1].set_xlabel('x [m]')
	ax[0].set_ylabel(' $\\xi$ [m]')
	ax[1].set_ylabel(' $\\nu$ [m/m] ')

	# place a text box in upper left in axes coords
	props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
	ax[0].text(0.75, 0.95, "$\omega =$" + str( omega )  , transform=ax[0].transAxes, fontsize=14, verticalalignment='top', bbox=props)
	pl.savefig('wave' + str(int(round(omega*10E45, 1))) + '.eps')

	omega *= 1.01

	pl.close('all')




# pl.figure()
# pl.plot(x,rho)
# pl.title('$\\rho$')
# pl.xlabel('x[m]')
# pl.ylabel('$\\rho$')
# pl.savefig('rho.eps')

# pl.figure()
# pl.plot(x,b*b)
# pl.title('$B^2$')
# pl.xlabel('x[m]')
# pl.ylabel('$B^2$')
# pl.savefig('b2.eps')

# pl.figure()
# pl.plot(x,vaSquared)
# pl.title('$v_a^2$')
# pl.xlabel('x[m]')
# pl.ylabel('$v_a^2$')
# pl.savefig('va2.eps')

# pl.figure()
# pl.plot(x,C1)
# pl.title('$C_1$')
# pl.xlabel('x[m]')
# pl.ylabel('$C_1$')
# pl.savefig('C1.eps')


# pl.show()
