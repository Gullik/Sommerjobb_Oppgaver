import sys
import numpy as np
import pylab as plt
from time import time
from mpl_toolkits.mplot3d import Axes3D

import sys
print (sys.version)

def dipoleField_rewritten(particle, nSteps, stepsize, method):
	
	e = 1.602E-19			#C
	Re = 6.371E6			#m
	mu0 = 4.*np.pi*1.E-7	#N/A^2

	steps = np.power(10, int(nSteps))
	timestep = np.power(10, int(stepsize))
	timestep = 1./timestep
	update = int(steps/10000)

	if particle == 'electron':
		q = - e
		m = 9.109E-31	#kg
	if particle == 'ion':
		q = e
		m = 16*1.674E-26  #kg

	C1 = q/m            #Used in the integration
	C2 = mu0/(4*np.pi)	#Constant in the magnetic field calculation

	#Setting initial conditions
	x = 5.*Re 		#m   
	y = 0.
	z = 0.

	vx = 300000.		#m/s
	vy = 0.
	vz = 300000.

	#Calculate suiting timesteps, let's assume the gyration frequency will stay at the same order as it is in the inital conditions
	Bx = bx(x,y,z,C2)
	By = by(x,y,z,C2)
	Bz = bz(x,y,z,C2)

	b_magnitude = r_length(Bx,By,Bz)
	# print b_magnitude
	period = np.abs(m/(q*b_magnitude)) * 2.*np.pi

	print("Period in s: " + str(period))

	#Open file
	results = open('results/' + particle + '_' + nSteps + '_'  + stepsize + '_' +   method +'.txt', 'w')

	print('Preparing to calculate the trajectory of an ' + particle + ' for ' + str(steps) + ' steps, steplength = ' +  str(timestep) + ' s, method ' + method)
	startTime = time()
	if method == 'Euler':

		i = 0
		while i < steps:
			#Magnetic field
			Bx = bx(x,y,z,C2)
			By = by(x,y,z,C2)
			Bz = bz(x,y,z,C2)

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

			if i%update == 0:

				#Calculating the parallel velocity by using dot(v,B)/|B| = v_parallel and the perpendicular velocity by using the cross(v,B)/|B| = v_perp
				b_magnitude = r_length(Bx, By, Bz)
				v_para = (vx*Bx + vy*By + vz*Bz)/b_magnitude
				v_perp = (vCrossB_x*vCrossB_x + vCrossB_y*vCrossB_y + vCrossB_z*vCrossB_z)
				v_perp = np.sqrt(v_perp)/b_magnitude
				kinetic_para = v_para*v_para
				kinetic_perp = v_perp*v_perp
				#Write to file
				results.write(str(x) + '\t' + str(y)  + ' \t'+ str(z) + '\t' + str(kinetic_para) + '\t' + str(kinetic_perp)  + '\t' + str(timestep*i) + '\n')
				if i % int(0.1*steps) == 0:
					print str( float(i) / float(steps) * 100)  + '%'

			i+=1

	elif method == 'Verlet':

		i = 0
		halfStep = timestep/2.

		while i < steps:
			#Magnetic field
			Bx = bx(x,y,z,C2)
			By = by(x,y,z,C2)
			Bz = bz(x,y,z,C2)

			#Update v(t + h/2)
			vx += halfStep*C1*crossx(vx,vy,vz,  Bx,By,Bz)
			vy += halfStep*C1*crossy(vx,vy,vz,  Bx,By,Bz)
			vz += halfStep*C1*crossz(vx,vy,vz,  Bx,By,Bz)

			#Update positions x(t+h)
			x += timestep*vx
			y += timestep*vy
			z += timestep*vz

			#Magnetic field
			Bx = bx(x,y,z,C2)
			By = by(x,y,z,C2)
			Bz = bz(x,y,z,C2)

			vCrossB_x = crossx(vx,vy,vz,  Bx,By,Bz)
			vCrossB_y = crossy(vx,vy,vz,  Bx,By,Bz)
			vCrossB_z = crossz(vx,vy,vz,  Bx,By,Bz)

			#Update v(t + h)
			vx += halfStep*C1*vCrossB_x 
			vy += halfStep*C1*vCrossB_y 
			vz += halfStep*C1*vCrossB_z


			if i%update == 0:
				#Calculating the parallel velocity by using dot(v,B)/|B| = v_parallel and the perpendicular velocity by using the cross(v,B)/|B| = v_perp
				b_magnitude = r_length(Bx, By, Bz)
				v_para = (vx*Bx + vy*By + vz*Bz)/b_magnitude
				v_perp = (vCrossB_x*vCrossB_x + vCrossB_y*vCrossB_y + vCrossB_z*vCrossB_z)
				v_perp = np.sqrt(v_perp)/b_magnitude
				kinetic_para = v_para*v_para
				kinetic_perp = v_perp*v_perp
				#Write to file
				results.write(str(x) + '\t' + str(y)  + ' \t'+ str(z) + '\t' + str(kinetic_para) + '\t' + str(kinetic_perp)  + '\t' + str(i*timestep) + '\n')
				if i % int(0.1*steps) == 0:
					print str( float(i) / float(steps) * 100)  + '%'

			i+=1

	endTime = time()
	print('Spent '+str(endTime - startTime) + ' s')

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
	# return 2.5*10**-7

def crossx(vx,vy,vz, Bx,By,Bz):
	#Returns the x component of the cross product between two vectors
	return vy*Bz - vz*By

def crossy(vx,vy,vz, Bx,By,Bz):	
	return vz*Bx - vx*Bz

def crossz(vx,vy,vz, Bx,By,Bz):
	return vx*By - vy*Bx

def analyze_data(name, method):
	print 'Analyzing data in ' + name

	data = np.genfromtxt('results/' + name + '.txt')

	x = data[:,0]
	y = data[:,1]
	z = data[:,2]
	E_para = data[:,3]
	E_perp = data[:,4]
	time= data[:,5]

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
	f.savefig('figures/' + name + 'xyz.eps')

	# print E_para
	f, ax = plt.subplots(3, sharex = True)
	plt.suptitle('Kinetic Energies, parallel and perpendicular to  $\\vec{B}$')
	ax[0].plot(time,E_para)
	ax[1].plot(time,E_perp)
	ax[2].plot(time, E_para + E_perp)
	ax[0].set_ylabel('$E_\parallel$')
	ax[1].set_ylabel('$E_\perp$')
	ax[2].set_ylabel('$E_{tot}$')
	ax[2].set_xlabel('Time [s]')

	f.savefig('figures/' + name + 'energy.eps')


	#3D plot of the trajectory
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	plt.suptitle('Particle trajectory')
	ax.plot(x,y,z)
	ax.set_xlabel('x [m]')
	ax.set_ylabel('y [m]')
	ax.set_zlabel('z [m]')

	ax.set_aspect('equal', 'datalim')

	fig.savefig('figures/' + name + '3Dplot.eps')

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
		

	plt.figure()
	plt.plot(time, phi*180/np.pi)

	print 'The drift period is: ' + str(  (2*np.pi) / (phi[-1] - phi[0]) * time[-1] )

	nBounce = 0
	for i in range(z.shape[0] - 1):
		if z[i] > 0 and z[i + 1] < 0:
			nBounce += 1

	print 'A rough estimate of the bounce period is: '	+ str( time[-1]/nBounce )


	plt.show()

if __name__ == '__main__':

	if sys.argv[5] == '1':
    		dipoleField_rewritten(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
		analyze_data(sys.argv[1] + '_' + sys.argv[2] + '_' + sys.argv[3] +'_' + sys.argv[4], sys.argv[4])
	else:
    		analyze_data(sys.argv[1] + '_' + sys.argv[2] + '_' + sys.argv[3] +'_' + sys.argv[4] , sys.argv[4])
    	




