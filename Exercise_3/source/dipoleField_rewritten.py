import sys
import numpy as np
import pylab as plt
from time import time
from mpl_toolkits.mplot3d import Axes3D

import sys
print (sys.version)

def dipoleField_rewritten(particle, steps):
	
	e = 1.602E-19			#C
	Re = 6.371E6			#m
	mu0 = 4.*np.pi*1.E-7	#N/A^2

	#Test new dipole
	# Re = 1.



	if particle == 'electron':
		q = - e
		m = 9.109E-31	#kg
	if particle == 'ion':
		q = e
		m = 16*1.674E-26  #kg

	steps = int(steps)

	C1 = q/m            #Used in the integration
	C2 = mu0/(4*np.pi)	#Constant in the magnetic field calculation

	#Setting initial conditions
	#	Since the changes in phasespace is quite small compared to the velocities and positions
	#	the phase space coordinates are shifted sligthly, r -> (-5Re, 0, 0) and v -> (0, 0, -300)
	# xShift = 5.*Re
	# vzShift = 300.

	x = 5.*Re 		#m   
	y = 0.
	z = 0.

	vx = 300.		#m/s
	vy = 0.
	vz = 1300.

	# dvx = 0.
	# dvy = 0.
	# dvz = 0.

	# v_para = 0
	# v_perp = 0

	# kinetic_para = 0
	# kinetic_perp = 0

	#Calculate suiting timesteps, let's assume the gyration frequency will stay at the same order as it is in the inital conditions
	Bx = bx(x,y,z,C2)
	By = by(x,y,z,C2)
	Bz = bz(x,y,z,C2)

	b_magnitude = r_length(Bx,By,Bz)

	period = np.abs(m/(q*b_magnitude)) * 2.*np.pi

	print("Period in s: " + str(period))

	timestep = period/100000.

	#Open file
	results = open('results_' + particle + '.txt', 'w')

	print('Preparing to calculate ' + str(steps) + ' steps, with an ' + particle)
	startTime = time()
	i = 0
	while i < steps:
		#Magnetic field
		Bx = bx(x,y,z,C2)
		By = by(x,y,z,C2)
		Bz = bz(x,y,z,C2)

		b_magnitude = r_length(Bx, By, Bz)
		vCrossB_x = crossx(vx,vy,vz,  Bx,By,Bz)
		vCrossB_y = crossy(vx,vy,vz,  Bx,By,Bz)
		vCrossB_z = crossz(vx,vy,vz,  Bx,By,Bz)


		#Update positions
		x += timestep*vx
		y += timestep*vy
		z += timestep*vz

		#Update v
		vx += timestep*C1*vCrossB_x 
		vy += timestep*C1*vCrossB_y 
		vz += timestep*C1*vCrossB_z

		# print i

		#Calculating the parallel velocity by using dot(v,B)/|B| = v_parallel and the perpendicular velocity by using the cross(v,B)/|B| = v_perp
		v_para = (vx*Bx + vy*By + vz*Bz)/b_magnitude
		v_perp = (vCrossB_x*vCrossB_x + vCrossB_y*vCrossB_y + vCrossB_z*vCrossB_z)
		v_perp = np.sqrt(v_perp)/b_magnitude
		kinetic_para = v_para*v_para
		kinetic_perp = v_perp*v_perp
		# print 'step ' + str(i)
		# print kinetic_para
		# print kinetic_perp
		# print vz

		# print np.sqrt(vx**2 + vy**2)

		if i%100000 == 0:
			# xx.append(x)
			# # yy.append(y)
			# zz.append(z)
			# # vvxx.append(vx)
			# # vvyy.append(vy)
			# # vvzz.append(vz)
			# K_para.append(v_para)
			# K_perp.append(v_perp)
			#Write to file
			results.write(str(x) + '\t' + str(y)  + ' \t'+ str(z) + '\t' + str(kinetic_para) + '\t' + str(kinetic_perp)  + '\n')
			if i % int(0.1*steps) == 0:
				print str( float(i) / float(steps) * 100)  + '%'

		i+=1

		

	endTime = time()
	print('Spent '+str(endTime - startTime) + ' s')

	# plot_trajectory(xx, zz, K_para, K_perp)

	return

def r_length(x,y,z):
	return (x*x + y*y +z*z)**.5

def dipMomentz(z):
	return 7.94E22

def bx(x,y,z, C2):
	mz = dipMomentz(z)
	return C2*(3.*x*(mz*z))/r_length(x,y,z)**5
	# return 0

def by(x,y,z,C2):
	mz = dipMomentz(z)
	return C2*(3.*y*(mz*z))/r_length(x,y,z)**5
	# return 0

def bz(x,y,z,C2):
	mz = dipMomentz(z)
	length = r_length(x,y,z)
	return C2*((3.*z*(mz*z))*length**-5 - mz*length**-3)


def crossx(vx,vy,vz, Bx,By,Bz):
	#Returns the x component of the cross product between two vectors
	return vy*Bz - vz*By

def crossy(vx,vy,vz, Bx,By,Bz):	
	return vz*Bx - vx*Bz

def crossz(vx,vy,vz, Bx,By,Bz):
	return vx*By - vy*Bx

def plot_from_file(name):
	data = np.genfromtxt(name)

	x = data[:,0]
	y = data[:,1]
	z = data[:,2]
	E_para = data[:,3]
	E_perp = data[:,4]
	timeArray = np.arange(x.shape[0])

	
	#Plotting positions
	f, ax = plt.subplots(3, sharex = True)
	
	ax[0].plot(timeArray,x)
	ax[1].plot(timeArray,y)
	ax[2].plot(timeArray,z)
	
	ax[0].set_ylabel('$x$')
	ax[1].set_ylabel('$y$')
	ax[2].set_ylabel('$z$')
	ax[2].set_xlabel('Timesteps')
	f.savefig('xyz.eps')


	f, ax = plt.subplots(3, sharex = True)
	plt.suptitle('Kinetic Energies, parallel and perpendicular to  $\\vec{B}$')
	ax[0].plot(timeArray,E_para)
	ax[1].plot(timeArray,E_perp)
	ax[2].plot(timeArray, E_para + E_perp)
	ax[0].set_ylabel('$E_\parallel$')
	ax[1].set_ylabel('$E_\perp$')
	ax[2].set_ylabel('$E_{tot}$')

	f.savefig('energy.eps')


	#3D plot of the trajectory
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	ax.plot(x,y, z)
	ax.set_xlabel('x')
	ax.set_ylabel('y')
	ax.set_zlabel('z')

	ax.set_aspect('equal', 'datalim')

	fig.savefig('3DPlot.eps')

	plt.show()





def plot_trajectory(xx, zz, K_para, K_perp):
	print( 'Reading data and preparing plots')

	startTime = time()

	# plt.figure()
	# plt.plot(xx, yy)
	# plt.title('Trajectory')
	# plt.xlabel('x [m]' )
	# plt.ylabel('y [m]')


	timeArray = range(0,len(xx))
	#Plotting positions
	f, ax = plt.subplots(2, sharex = True)
	
	ax[0].plot(timeArray,xx)
	ax[1].plot(timeArray,yy)
	ax[2].plot(timeArray,zz)
	
	ax[0].set_ylabel('$x$')
	ax[1].set_ylabel('$y$')
	ax[2].set_ylabel('$z$')
	ax[2].set_xlabel('Timesteps')
	f.savefig('xyz.eps')

	# #Plotting velocities
	# f, ax = plt.subplots(3, sharex = True)
	
	# ax[0].plot(timeArray,vvxx)
	# ax[1].plot(timeArray,vvyy)
	# ax[2].plot(timeArray,vvzz)
	
	# ax[0].set_ylabel('$v_x$')
	# ax[1].set_ylabel('$v_y$')
	# ax[2].set_ylabel('$v_z$')
	# ax[2].set_xlabel('Timesteps')

	# f.savefig('v_xyz.eps')

	f, ax = plt.subplots(2, sharex = True)
	plt.suptitle('Kinetic Energies, parallel and perpendicular to  $\\vec{B}$')
	ax[0].plot(timeArray,K_para)
	ax[1].plot(timeArray,K_perp)

	

	f.savefig('energy.eps')

	# #3D plot of the trajectory
	# fig = plt.figure()
	# ax = fig.add_subplot(111, projection='3d')

	# ax.plot(xx,yy, zz)
	# ax.set_xlabel('x')
	# ax.set_ylabel('y')
	# ax.set_zlabel('z')

	# ax.set_aspect('equal')

	# fig.savefig('3DPlot.eps')


	endTime = time()

	print ('Spent ' + str(endTime - startTime) + ' s')
	plt.show()
	


if __name__ == '__main__':

	if sys.argv[3] == '1':
    		dipoleField_rewritten(sys.argv[1], sys.argv[2])
		plot_from_file('results_' + sys.argv[1] +'.txt')
	else:
    		plot_from_file('results_' + sys.argv[1] +'.txt')
    	




