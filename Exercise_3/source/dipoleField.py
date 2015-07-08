import numpy as np
import scipy as sp
import pylab as pl
import scipy.stats as stats
from mpl_toolkits.mplot3d import Axes3D


#Declaration of several variables, all the units are in SI
e = -1.602E-19		#C
m_e = 9.109E-31		#kg
m_i = 16*1.674E-26  #kg
Re = 6.371E3		#m
mu = 4.*np.pi*1.E-7	#N/A^2
dipoleMoment = 7.94 *1.E22		#A m^2	Dipole moment

particle = "electron"
# particle = "ion"

# #Choose if for ion or electron
if particle == "electron":
	# #electron settings
	h = 1.E-22		#stepsize
	m = m_e	
	q = e

if particle == "ion":
	#Ion settings
	h = 1.E-8
	m = m_i
	q = -2*e

steps = int(1E8)

#Defining the vectors which contains the positions and velocities during the simulation
r = np.zeros((steps,3))	#x, y and z direction
# yy = np.zeros(steps)
# zz = np.zeros(steps)
v = np.zeros((steps,3))	#velocities in x y and z directions
# vv = np.zeros(steps)
# ww = np.zeros(steps)
kineticEnergy = np.zeros(steps)
kineticEnergyParalell = np.zeros(steps)
kineticEnergyPerpendicular = np.zeros(steps)

#InitialConditions
X = 5*Re	#m
Y = 0.
Z = 0.
U = 300.	#m/s speed in x direction
V = 0.		#m/s velocity in y direction	
W = 300.	#m/S velocity in z direction

r[0,:] = [X,Y,Z]
v[0,:] = [U,V,W]


C2 = 3*mu/(4*np.pi)
def magneticField(r):
	distance = np.sqrt(r[0]*r[0]+ r[1]*r[1] + r[2]*r[2])
	return (C2*r*(r[2]*m)/distance**5) - [0,0,dipoleMoment]/(distance**3)
	# return [0,0,dipoleMoment]/distance**3


# print np.cross( v[0,:],magneticField(r[0,:]) )
# print q*np.linalg.norm(magneticField(r[0,:]))/m
# print np.linalg.norm(magneticField(r[0,:]))
# distance = np.sqrt(X*X + Y*Y + Z*Z)
# print dipoleMoment
# print distance
# print distance**3
# print dipoleMoment/distance**3
# print magneticField(r[0,:])
# print 'Gyrofrequency: ' + str(q*np.linalg.norm(magneticField(r[0,:]))/m)

#Precalculating some constants so not to calculate them in the loop
C1 = h*(q/m)

print C1

for i in range(0,steps -1):
	#Rewrote it so it doesn't go searching through the vectors so often, faster

	#Let us first calculate the steps since last move
	B = magneticField(r[i,:])
	# print 'New step'
	# print B
	dv = C1 * np.cross(v[i,:],B)
	# print dv
	# print v[i,:]
	# dv = C1 * np.cross(v[i,:],B0)

	dr = h*v[i,:]

	v[i + 1,:] = v[i,:] + dv
	r[i + 1,:] = r[i,:] + dr

	# #Calculate new positions and velocities
	# X = X + dX
	# Y = Y + dY
	# U = U + dU
	# # V = V + dV
	# W = W +dW

	#Storing datapoints
	# uu[i + 1] = U
	# vv[i + 1] = V
	# ww[i + 1] = W
	# xx[i + 1] = X
	# yy[i + 1] = Y
	# zz[i + 1] = Z 
	if i % int(0.25*steps) == 0:
		print str( float(i) / float(steps) * 100)  + '%'

# print v

fig = pl.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot(r[:,0],r[:,1], r[:,2])
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
# ax.set_aspect('equal', 'datalim')


# #Trying to plot the velocities
# time = range(0,steps)
# pl.figure()
# pl.plot(time, v[:,0])
# pl.figure()
# pl.plot(time, v[:,1])
# pl.figure()
# pl.plot(time, v[:,2])


pl.show()


