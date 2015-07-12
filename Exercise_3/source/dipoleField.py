import numpy as np
import scipy as sp
import pylab as pl
import scipy.stats as stats
from mpl_toolkits.mplot3d import Axes3D
from time import time


#Declaration of several variables, all the units are in SI
e = -1.602E-19		#C
m_e = 9.109E-31		#kg
m_i = 16*1.674E-26  #kg
Re = 6.371E3		#m
mu = 4.*np.pi*1.E-7	#N/A^2
dipoleMoment = 7.94 *1.E22		#A m^2	Dipole moment  	Note, this is directed along the z-axis
dipoleVector = np.zeros(3)
dipoleVector[2] = dipoleMoment


# particle = "electron"
particle = "ion"

# #Choose if for ion or electron
if particle == "electron":
	# #electron settings
	h = 1.E-23		#stepsize
	m = m_e	
	q = e

if particle == "ion":
	#Ion settings
	h = 1.E-19
	m = m_i
	q = -2*e

steps = int(1E6)

#Defining the vectors which contains the positions and velocities during the simulation
r = np.zeros((steps,3))	#x, y and z direction
v = np.zeros((steps,3))	#velocities in x y and z directions


kineticEnergy = np.zeros(steps)
kineticEnergyParalell = np.zeros(steps)
kineticEnergyPerpendicular = np.zeros(steps)
v_para = np.zeros(steps)
v_perp = np.zeros(steps)
v_magn = np.zeros(steps)

#InitialConditions 	#Note! The coordinate system was shifted 5*Re and the velocity was shifted 300 in v_z axis
X = 0 #5*Re	#m
Y = 0.
Z = 0.
U = 300.	#m/s speed in x direction
V = 0.		#m/s velocity in y direction	
W = 0.		#m/S velocity in z direction

xShift = 5.*Re 			#Shift in the x-axis to avoid rounding away the results
vzShift = 3.E2			#Shift in the velocity axis to avoid rounding away the results

r[0,:] = [X,Y,Z]	
v[0,:] = [U,V,W]
R = np.zeros((1,3))	#R and V is used in the for loop to avoid going into the vectors more than necessary
V = np.zeros((1,3))
r2 = np.zeros(3)	#Used to compensate for the shifted axis in the magnetic field
MagField = np.zeros((steps,3)) #Storing the field values


C1 = h*(q/m)			#Constant needed in the calculation of the new velocity
C2 = 3.*mu/(4.*np.pi)	#Constant used in the magnetic field calculation


def magneticField(r):
	#Moved the coordinate system to insure that the movement in the x-axis is not rounded away

	r2[0] = r[0] + xShift
	r2[1] = r[1]
	r2[2] = r[2]

	distance = np.sqrt(r2[0]*r2[0] + r2[1]*r2[1] + r2[2]*r2[2])

	return (C2*r2*(r2[2]*dipoleMoment)/distance**5) - dipoleVector/(distance**3)

start = time()
for i in range(0,steps -1):
	R = r[i,:]
	V = v[i,:]

	B = magneticField(R)
	dv = C1 * np.cross(V ,B)
	dr = h*(V + [0,0,vzShift])

	print dr[2]/vzShift

	v[i + 1,:] = V + dv
	r[i + 1,:] = R + dr

	#The magnetic field is also stored since we need it later when calculating the paralell and perpendilar velocity
	MagField[i,:] = B


	if i % int(0.25*steps) == 0:
		print str( float(i) / float(steps) * 100)  + '%'
end = time()

print 'Time used: ' + str(end - start)


fig = pl.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot(r[:,0],r[:,1], r[:,2])
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

ax.set_aspect('equal')

# #Trying to plot the velocities
time = range(0,steps)
f, axArray = pl.subplots(3, sharex = True)

d = 0
for ax in axArray.flat:
	ax.plot(time,v[:,d])
	d = d+1

axArray[0].set_ylabel('$v_x$')
axArray[1].set_ylabel('$v_y$')
axArray[2].set_ylabel('$v_z$')

#Calculate the perpendicular and parallel velocity
#Normalize the magnetic field so we can calculate v_para = v . b



MagField[steps-1,:] = magneticField(r[steps-1,:])	#Adding last entry, wasn't done in the for-loop
norm = np.linalg.norm(MagField, axis = -1)
MagField/=norm[:,None]

v_magn = np.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])

for i in range(0,steps ):
	v_para[i] = np.dot(v[i], MagField[i])

# v_perp = v_magn[:] - v_para

# print v_perp
# print v_para

f, axArray = pl.subplots(1, sharex = True)

# axArray[0].plot(time,v_perp)
axArray.plot(time,v_para)

# print v_magn
# print v_para




# #Plotting the positions

# f, axArray = pl.subplots(3, sharex = True)

# d = 0
# for ax in axArray.flat:
# 	ax.plot(time,r[:,d])
# 	d = d+1
# axArray[0].set_ylabel('$x$')
# axArray[1].set_ylabel('$y$')
# axArray[2].set_ylabel('$z$')

# pl.figure()
# pl.plot(time ,r[:,1])

pl.show()


