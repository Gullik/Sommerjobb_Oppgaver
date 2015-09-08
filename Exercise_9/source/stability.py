import numpy as np
import pylab as plt
from mpl_toolkits.mplot3d import Axes3D
import sys as sys
from time import time


# This function calculates the trajectory of a particle, being affected by the earths and gravity
# It is initialized by for example
# python stability.py 5 3 Euler



def r_length(x,y,z):
	return np.sqrt(x*x + y*y +z*z)

def dipMomentz():
	return 7.94E22			#A m^2

def bx(x,y,z, C2):
	mz = dipMomentz()
	return C2*(3.*x*(mz*z))/r_length(x,y,z)**5

def by(x,y,z,C2):
	mz = dipMomentz()
	return C2*(3.*y*(mz*z))/r_length(x,y,z)**5

def bz(x,y,z,C2):
	mz = dipMomentz()
	length = r_length(x,y,z)
	return C2*((3.*z*z*mz)*length**-5 - mz*length**-3)

def b_magnitude(x,y,z, C):

	Bx = bx(x,y,z, C)
	By = by(x,y,z, C)
	Bz = bz(x,y,z, C)

	return r_length(Bx,By,Bz)

def gx(x,y,z, C):
	return C*x*r_length(x,y,z)**-3

def gy(x,y,z, C):
	return C*y*r_length(x,y,z)**-3

def gz(x,y,z, C):
	return C*z*r_length(x,y,z)**-3

def dB_zdz(z, C):
	return -C*z**-4

def crossx(vx,vy,vz, Bx,By,Bz):
	#Returns the x component of the cross product between two vectors
	return vy*Bz - vz*By

def crossy(vx,vy,vz, Bx,By,Bz):	
	return vz*Bx - vx*Bz

def crossz(vx,vy,vz, Bx,By,Bz):
	return vx*By - vy*Bx

def b(r,C2):
	#Takes all the magnetic field and returns them as a 3d vector
	x = r[0]
	y = r[1]
	z = r[2]

	return np.array([bx(x,y,z,C2), by(x,y,z,C2), bz(x,y,z,C2)])

def gravity(r, C):
	#Takes all the gravitational force and returns them as a 3d vector
	x = r[0]
	y = r[1]
	z = r[2]

	return np.array([gx(x,y,z,C), gy(x,y,z,C), gz(x,y,z,C)])


def force(x,y,z, vx,vy,vz, C1, C2, C3):
	#Calculates the magnetic force and the gravitational force for a particle
	Bx = bx( x,y,z , C2)
	By = by( x,y,z , C2)
	Bz = bz( x,y,z , C2)

	force_x = crossx(vx,vy,vz, Bx, By,Bz)*C1
	force_y = crossy(vx,vy,vz, Bx, By,Bz)*C1
	force_z = crossz(vx,vy,vz, Bx, By,Bz)*C1
	
	force_x -= gx( x,y,z , C3)
	force_y -= gy( x,y,z , C3)
	force_z -= gz( x,y,z , C3)

	# print  force_z

	return force_x, force_y, force_z

def forceRK4( y, mag_field, C1, C2, C3 ):
	#	This function takes as the input a phase vector of positions, and velocities, and a magnetic field.
	# 	Then it calculates the the function f = (v, q/m*vxB) and returns that phase vector
	temp = np.zeros(6)

	temp[:3] = y[3:]										#Calculation of positions dr/dt = v
	temp[3:] = C1*np.cross( y[3:], mag_field(y[:3], C2) )	#Calculation of force 	  dv/dt	= q/m*(v x B)
	temp[3:] -= gravity(y[:3], C3) 

	return temp

def gyration_radius(v_perp, B, m, q):
	#Returns the gyration radius rho = m*v_perp/(qB)
	return np.abs(m*v_perp/(q*B))

def rk4(y, dydt, C1, C2, C3, steps, timestep):
	#This takes a function and it's derivative and
	#uses the RK4 method, and returns a vector of positions

	# C1 = q/m
	# C2 = mu_0/(4pi)
	# C3 = gamma*M_earth
 
	k1 = np.zeros(6)
	k2 = np.zeros(6)
	k3 = np.zeros(6)
	k4 = np.zeros(6)

	update = 1
	# update = int(steps/1000)
	data = np.zeros((6,int( steps/update)))

	for i in range(steps):

		print "Run " + str(i)

		print forceRK4(y, b, C1,C2,C3)[5]


		k1 = timestep*dydt( y       , b, C1, C2, C3)
		k2 = timestep*dydt( y + k1/2, b, C1, C2, C3)
		k3 = timestep*dydt( y + k2/2, b, C1, C2, C3)
		k4 = timestep*dydt( y + k3  , b, C1, C2, C3)

		print 'dz =' + str( (k1 + 2*k2 + 2*k3 + k4 )[2]/6.) +  \
				'\t dv_z =' + str( (k1 + 2*k2 + 2*k3 + k4 )[5]/6.)

		y += ( k1 + 2*k2 + 2*k3 + k4 )/6.
		# y += k1

		#Store the present data
		# # if i%update == 0:
		# 	data[:,int(i/update)] = y
		data[:,int(i/update)] = y

	return data

def euler(force, x, y, z, vx, vy, vz,C1, C2, C3, nSteps, timestep ):

	# C1 = q/m
	# C2 = mu_0/(4pi)
	# C3 = gamma*M_earth
	
	positions = np.zeros((nSteps,3))
	velocities = np.zeros((nSteps, 3))

	i = 0
	while i < nSteps:

		force_x, force_y, force_z = force( x,y,z, vx,vy,vz, C1, C2, C3)

		#Update positions
		x += timestep*vx
		y += timestep*vy
		z += timestep*vz

		#Update v
		vx += timestep*force_x
		vy += timestep*force_y 
		vz += timestep*force_z

		positions[i,:] = [x,y,z]
		velocities[i,:] = [vx,vy,vz]

		i += 1

	#Plots
	time = np.linspace(0, nSteps*timestep, nSteps)

	return time, positions, velocities

def plot_routine( time, positions, velocities ):

	# #3D plot of the trajectory
	# fig = plt.figure()
	# ax = fig.add_subplot(111, projection='3d')
	# plt.suptitle('Particle trajectory')
	# ax.plot(positions[:,0],positions[:,1],positions[:,2])
	# ax.set_xlabel('x [m]')
	# ax.set_ylabel('y [m]')
	# ax.set_zlabel('z [m]')

	# ax.set_aspect('equal', 'datalim')

	# fig.savefig('figures/3Dplot.eps')

	fig = plt.figure()
	ax = fig.add_subplot(111)
	plt.plot(time,velocities[:,2])
	plt.xlabel('Time $[s]$')
	plt.ylabel('Vertical Velocity $[m/s]$')

	fig.savefig('figures/vertical_vel.eps')

	fig = plt.figure()
	ax = fig.add_subplot(111)
	plt.plot(time,positions[:,2])
	plt.xlabel('Time $[s]$')
	plt.ylabel('z $[m/s]$')

	print np.max(positions[:,2]) - np.min(positions[:,2])

	fig.savefig('figures/vertical_pos.eps')

	#Plotting positions
	f, ax = plt.subplots(3, sharex = True)
	plt.suptitle('Positions')

	
	ax[0].plot(time,positions[:,0])
	ax[1].plot(time,positions[:,1])
	ax[2].plot(time,positions[:,2])
	
	ax[0].set_ylabel('$x [m]$')
	ax[1].set_ylabel('$y [m]$')
	ax[2].set_ylabel('$z [m]$')
	ax[2].set_xlabel('Time [s]')
	f.savefig('figures/xyz.eps')

	return 0 #positions

def initial_velocity(z, pitch, C1, C2, C3, C4):

	# #Calculating a velocity that will balance the forces, along with initialplacement
	# tol = 1.E-10
	x = 0
	y = 0

	#Calculating a velocity that will balance the forces
	v_perp = np.abs(2*b_magnitude(0,0,z, C2)*gz(0,0,z,C3)/dB_zdz(z,C4))

	v_perp = np.sqrt(v_perp)

	x = gyration_radius(v_perp, b_magnitude(x,y,z, C2), m_e, q )
	y = 0
	
	vx = v_perp*np.cos(pitch)
	vy = v_perp*np.sin(pitch)
	vz = int(sys.argv[2])

	# force_x, force_y, force_z = force( x,y,z, vx + 10000,vy,vz, C1, C2, C3)
	print forceRK4([x,y,z,vx,vy,vz], b, C1,C2,C3)[5]
	
	return x, y, z, vx, vy, vz

if __name__ == '__main__':

	mu0 	= 4.*np.pi*1.E-7	#N/A^2
	Re 		= 6.371E6			#m
	q 		= -1.602E-19		#C
	m_e 	= 9.109E-31			#kg	
	# m_i 	= 
	Me 		= 5.9E24			#kg
	gamma 	= 6.674E-11			#Nm^2/kg^2 
	pitch	= np.deg2rad(90)	#radian

	nSteps 	= np.power(10,int(sys.argv[1]))

	if len(sys.argv) == 1:
		print "The arguments for steps, initial vertical velocity and method is missing"
		print " python stability.py steps v_z Euler/rk4 "
	
	C1 = q/m_e
	C2 = mu0/(4*np.pi)	#Constant in the magnetic field calculation
	C3 = gamma*Me
	C4 = 3*mu0*dipMomentz()/(2*np.pi)
	

	#Initial position
	z = 200.E3 + Re				#m

	x, y, z, vx, vy, vz = initial_velocity(z, pitch, C1, C2, C3, C4)

	timeStart = time()
	if sys.argv[3] == "Euler":

		timestep 	= 1.E-10
		t, positions, velocity = euler(force, x, y, z, vx, vy, vz, C1, C2, C3, nSteps, timestep )

	elif sys.argv[3] == "rk4":

		timestep 	= 1.E-8

		dydt = forceRK4
		y = np.array( (x, y, z, vx, vy, vz) )

		data = rk4(y, dydt, C1, C2, C3, nSteps, timestep)

		positions = data[:3,:]
		velocity = data[3:,:]

		positions = positions.transpose()
		velocity = velocity.transpose()

		t = np.arange(positions.shape[0])*timestep
	else:
		print "No method chosen"

	timeEnd = time()

	print "The simulation ran for "   + str(timestep*nSteps) +     \
			"s, and used " + str( timeEnd - timeStart ) + 's'

	plot_routine(t, positions, velocity)

	plt.show()

