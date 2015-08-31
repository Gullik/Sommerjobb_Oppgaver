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

	# print np.array([bx(x,y,z,C2), by(x,y,z,C2), bz(x,y,z,C2)])
	 
	return np.array([bx(x,y,z,C2), by(x,y,z,C2), bz(x,y,z,C2)])


def function(y, mag_field, C1, C2):
	#	This function takes as the input a phase vector of positions, and velocities, and a magnetic field.
	# 	Then it calculates the the function f = (v, q/m*vxB) and returns that phase vector
	temp = np.zeros(6)

	temp[:3] = y[3:]														#Calculation of positions dr/dt = v
	temp[3:] = C1*np.cross( y[3:], mag_field(y[:3], C2) )					#Calculation of force 	  dv/dt	= q/m*(v x B)
	
	return temp


def rk4(y, dydt, C1, C2, steps, timestep):
	#This takes a function and it's derivative and
	#uses the RK4 method, and returns a vector of positions

	k1 = np.zeros(6)
	k2 = np.zeros(6)
	k3 = np.zeros(6)
	k4 = np.zeros(6)

	# update = 250
	update = int(steps/10000)
	data = np.zeros((8,int( steps/update)))

	for i in range(steps):

		k1 = timestep*dydt( y       , b, C1, C2 )
		k2 = timestep*dydt( y + k1/2, b, C1, C2 )
		k3 = timestep*dydt( y + k2/2, b, C1, C2 )
		k4 = timestep*dydt( y + k3  , b, C1, C2 )

		y += ( k1 + 2*k2 + 2*k3 + k4 )/6.

		#Store the present data

		if i%update == 0:
			data[:-2,int(i/update)] = y

			#Calculating the perpendicular and parallel energy
			b_field = b(y[:3], C2)
			b_length = r_length(b_field[0], b_field[1], b_field[2])

			v_para = (y[3]*b_field[0] +y[4]*b_field[1] + y[5]*b_field[2])/b_length
			v_perp = np.cross(y[3:],b_field)/b_length
			
			data[6,int(i/update)] = v_para*v_para
			data[7,int(i/update)] = np.dot(v_perp,v_perp)


	return data

def plots(data, timestep, mass):

	x = data[0,:]
	y = data[1,:]
	z = data[2,:]

	vx = data[3,:]
	vy = data[4,:]
	vz = data[5,:]

	E_para = mass/2.*data[6,:]
	E_perp = mass/2.*data[7,:]
	E_tot  = mass/2.*(vx*vx + vy*vy + vz*vz)

	#Converting energy to eV
	eV = 6.24E18
	E_para *= eV
	E_perp *= eV
	E_tot  *= eV

	time = np.arange(data.shape[1])*timestep

	#Plotting positions
	f, ax = plt.subplots(3, sharex = True)
	plt.suptitle('Positions')

	ax[0].plot(time,x)
	ax[1].plot(time,y)
	ax[2].plot(time,z)
	
	ax[0].set_ylabel('$x [m]$')
	ax[1].set_ylabel('$y [m]$')
	ax[2].set_ylabel('$z [m]$')
	ax[2].set_xlabel('Time [s]')

	f.savefig('figures/rk4_xyz.eps')

	# f, ax = plt.subplots(3, sharex = True)
	# plt.suptitle('Velocitites')

	# ax[0].plot(time,vx)
	# ax[1].plot(time,vy)
	# ax[2].plot(time,vz)
	
	# ax[0].set_ylabel('$v_x [m/s]$')
	# ax[1].set_ylabel('$v_y [m/s]$')
	# ax[2].set_ylabel('$v_z [m/s]$')
	# ax[2].set_xlabel('Time [s]')

	# f.savefig('figures/rk4_v.eps')

	#3D plot of the trajectory
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	plt.suptitle('Particle trajectory')
	ax.plot(x,y,z, lw = 0.4, color = 'r')
	ax.set_xlabel('x [m]')
	ax.set_ylabel('y [m]')
	ax.set_zlabel('z [m]')

	ax.set_aspect('equal', 'datalim')

	fig.savefig('figures/rk4_3D.eps')

	#Plotting the energy, 1/2 mv^2
	f, ax = plt.subplots(3, sharex = True)
	plt.suptitle('Energy')

	ax[0].plot(time,E_tot)
	ax[1].plot(time,E_para)
	ax[2].plot(time,E_perp)
	
	ax[0].set_ylabel('$E [eV]$')
	ax[1].set_ylabel('$E_{\parallel} [eV]$')
	ax[2].set_ylabel('$E_{\perp} [eV]$')
	ax[2].set_xlabel('Time [s]')

	f.savefig('figures/rk4_E.eps')


	fig = plt.figure()
	ax = fig.add_subplot(111)
	plt.suptitle('Total Energy $[eV]$')
	ax.plot(time,E_tot)
	ax.set_xlabel('Time [s]')
	ax.set_ylabel('$E$')

	fig.savefig('figures/rk4Energy.eps')

	# #Calculating the gyro period, bounce and drift period
	phi = np.arctan2(x,y)

	#Calculating gyrations as fast periodic changes in the xy angle, only works if the changes in angle is very stable
	nGyro = 0
	up = False
	for i in range(phi.shape[0] - 1):
		if up:
			if phi[i] > phi[i + 1]:
				nGyro += 1
				up = False
		else:
			if phi[i] < phi[i + 1]:
				up = True

	print 'The average gyration period is: ' + str(time[-1]/nGyro)
		
	print 'The drift period is: ' + str(  (2*np.pi) / (phi[-1] - phi[0]) * time[-1] )

	nBounce = 0
	for i in range(z.shape[0] - 1):
		if z[i] > 0 and z[i + 1] < 0:
			nBounce += 1

	print 'A rough estimate of the bounce period is: '	+ str( time[-1]/nBounce )

	return



if __name__ == '__main__':

	#Defining the physical parameters
	e = 1.602E-19			#C
	Re = 6.371E6			#m
	mu0 = 4.*np.pi*1.E-7	#N/A^2
	q = e
	m = 16*1.674E-26  		#kg

	steps = 3*10**int(sys.argv[1])
	timestep = 10**int(sys.argv[2])

	# steps = 1

	C1 = q/m            #Used in the integration
	C2 = mu0/(4*np.pi)	#Constant in the magnetic field calculation
	
	#Setting initial conditions
	r = np.array([5.*Re, 0., 0.]) 		#m   

	pitch 	= np.deg2rad(25)
	v_0 	= 10**5

	v = np.array([v_0*np.cos(pitch), 0.,  v_0*np.sin(pitch)])		#m/s

	y = np.zeros(6)
	y = np.array([r[0], r[1], r[2], v[0], v[1], v[2]])

	dydt = function

	startTime = time()
	data = rk4(y, function, C1, C2, steps, timestep)
	endTime = time()

	print "Used " + str( endTime - startTime ) + 's'

	plots(data, timestep, m)

	plt.show()

