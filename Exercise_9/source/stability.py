import numpy as np
import pylab as plt
from mpl_toolkits.mplot3d import Axes3D
import sys as sys



def r_length(x,y,z):
	return np.sqrt(x*x + y*y +z*z)

def dipMomentz():
	return 7.94E22			#A m^2

def bx(x,y,z, C2):
	mz = dipMomentz()
	return C2*(3.*x*(mz*z))/r_length(x,y,z)**5
	# return 0

def by(x,y,z,C2):
	mz = dipMomentz()
	return C2*(3.*y*(mz*z))/r_length(x,y,z)**5
	# return 0

def bz(x,y,z,C2):
	mz = dipMomentz()
	length = r_length(x,y,z)
	return C2*((3.*z*z*mz)*length**-5 - mz*length**-3)

def b_magnitude(x,y,z, C):

	Bx = bx(x,y,z, C)
	By = by(x,y,z, C)
	Bz = bz(x,y,z, C)

	# print Bz

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

def force(x,y,z, vx,vy,vz, C1, C2, C3, C4):
	#Calculates the magnetic force and the gravitational force for a particle
	Bx = bx( x,y,z , C2)
	By = by( x,y,z , C2)
	Bz = bz( x,y,z , C2)

	# print Bx
	# print By

	# force_z = 0


	force_x = crossx(vx,vy,vz, Bx, By,Bz)*C4
	force_y = crossy(vx,vy,vz, Bx, By,Bz)*C4
	force_z = crossz(vx,vy,vz, Bx, By,Bz)*C4
	
	# print force_z

	# print r_length(vx,vy,vz)
	# print r_length(Bx,By,Bz)
	# print 'Magnetic force:', force_z
	# print 'Graviational force: ', gz( x,y,z , C1)
	# print 'ratio of the forces', force_z/gz( x,y,z , C1)

	force_x -= gx( x,y,z , C1)
	force_y -= gy( x,y,z , C1)
	force_z -= gz( x,y,z , C1)

	# print force_z

	return force_x, force_y, force_z

def gyration_radius(v_perp, B, m, q):
	#Returns the gyration radius rho = m*v_perp/(qB)
	return np.abs(m*v_perp/(q*B))

def euler(force, x, y, z, vx, vy, vz,C1, C2, C3, C4, nSteps, timestep ):
	
	positions = np.zeros((nSteps,3))
	vert_force = np.ones(nSteps)


	i = 0
	while i < nSteps:


		force_x, force_y, force_z = force( x,y,z, vx,vy,vz, C1, C2, C3, C4)

		#Update positions
		x += timestep*vx
		y += timestep*vy
		z += timestep*vz

		#Update v
		vx += timestep*force_x
		vy += timestep*force_y 
		vz += timestep*force_z

		positions[i,:] = [x,y,z]
		vert_force[i] = force_z

		i += 1

	#Plots
	time = np.linspace(0, nSteps*timestep, nSteps)

	# print positions[0,2]
	# print positions[-1,2]
	
	#3D plot of the trajectory
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	plt.suptitle('Particle trajectory')
	ax.plot(positions[:,0],positions[:,1],positions[:,2])
	ax.set_xlabel('x [m]')
	ax.set_ylabel('y [m]')
	ax.set_zlabel('z [m]')

	# ax.set_aspect('equal', 'datalim')

	fig.savefig('../figures/3Dplot.eps')


	fig = plt.figure()
	ax = fig.add_subplot(111)
	plt.plot(time,vert_force)


	return 0 #positions

if __name__ == '__main__':

	mu0 	= 4.*np.pi*1.E-7	#N/A^2
	Re 		= 6.371E6			#m
	q 		= -1.602E-19		#C
	m_e 	= 9.109E-31			#kg	
	Me 		= 5.9E24			#kg
	gamma 	= 6.674E-11			#Nm^2/kg^2 
	pitch	= np.deg2rad(90)	#radian

	timestep 	= 1.E-12
	nSteps 		= np.power(10,int(sys.argv[1]))

	C1 = gamma*Me
	C2 = mu0/(4*np.pi)	#Constant in the magnetic field calculation
	C3 = 3*mu0*dipMomentz()/(2*np.pi)
	C4 = q/m_e

	#Initial position
	z = 200.E3 + Re				#m

	#Calculating a velocity that will balance the forces
	v_perp = np.abs(2*b_magnitude(0,0,z, C2)*gz(0,0,z,C1)/dB_zdz(z,C3))

	v_perp = np.sqrt(v_perp)

	# print 'g Field: ' + str(gz(0,0,z,C1))
	# print 'B Field: ' + str(b_magnitude(0,0,z, C2))
	# print dB_zdz(z,C3)
	# print C3


	# print m_e*v_perp*dB_zdz(z,C3)/(2*b_magnitude(0,0,z, C2))
	# print m_e*gz(0,0,z,C1)

	x = gyration_radius(v_perp, b_magnitude(0,0,z, C2), m_e, q )
	y = 0

	print v_perp

	vx = v_perp*np.cos(pitch)
	vy = v_perp*np.sin(pitch)
	vz = int(sys.argv[2])

	euler(force, x, y, z, vx, vy, vz, C1, C2, C3, C4, nSteps, timestep )

	plt.show()

