import numpy as np
import scipy as sp
import pylab as pl
import scipy.stats as stats



#Declaration of several variables, all the units are in SI
e =  1.602E-19		#C
m_e = 9.109E-31		#kg
m_i = 16*1.674E-27 	#kg
B0 =  50000E-9		#Tesla


#Choose if for ion or electron
#electron settings
h = 1.E-12
m = m_e
q = -e

# #Ion settings
# h = 1.*1.E-7
# m = m_i
# q = e

steps = int(1E7)

#Defining the vectors which contains the positions and velocities during the simulation
x = np.zeros(steps)
y = np.zeros(steps)
vx = np.zeros(steps)
vy = np.zeros(steps)
vx_half = np.zeros(steps)
vy_half = np.zeros(steps)
kineticEnergy = np.zeros(steps)
angles = np.zeros(steps)
time = np.linspace(0,steps*h, steps) #Only used for plotting the kinetic energy

#Setting initial conditions
x[0] = 0
y[0] = 0
vx[0] = 500
vy[0] = 0

#Precalculating some constants so not to calculate them in the loop
C1 = h*B0*q/m

# print C1

for i in range(0,steps -1):
	x[i + 1] = x[i] + h*vx[i]
	y[i + 1] = y[i] + h*vy[i]
	vx[i + 1] = vx[i] + C1*vy[i]
	vy[i + 1] = vy[i] - C1*vx[i]

	# vx_half[i +1] = vx[i] + C1/2*vy[i]
	# vy_half[i + 1] = vy[i] - C1/2*vx[i]

	# x[i + 1] = x[i] + h*vx_half[i]
	# y[i + 1] = y[i] + h*vy_half[i]
	# vx[i + 1] = vx_half[i] + C1/2*vy_half[i]
	# vy[i + 1] = vy_half[i] - C1/2*vx_half[i]


	if i % int(0.25*steps) == 0:
		print str( float(i) / float(steps) * 100)  + '%'


#Calculating the kinetic energy to see that it stays decently constant
kineticEnergy = m*(vx*vx + vy*vy)/2.
print "relative change E_K = " + str((np.max(kineticEnergy) - np.min(kineticEnergy))/np.mean(kineticEnergy)) 
print "std(E_K) = " + str(stats.tstd(kineticEnergy))

pl.figure()
pl.plot(x,y)
pl.title('Gyration of an electron')
# pl.title('Gyration of an oxygen ion')
pl.xlabel("x [m]")
pl.ylabel("y [m]")
pl.axes().set_aspect('equal', 'datalim')

# pl.figure()
# time = np.arange(kineticEnergy.shape[0])
# pl.plot(time,kineticEnergy)

pl.savefig("electronGyration.eps")
# pl.savefig("ionGyration.eps")




# pl.figure()
# pl.plot(time, kineticEnergy)		#This wasn't a really good plot for anything

# pl.show()

#Time to calculate the gyration period and gyration radius
center = (np.mean(x), np.mean(y))

rho = np.mean( np.sqrt((x-center[0])**2 + (y-center[1])**2 ) )

rotations = 0
lastRotationTime = 0
startRotationTime = 0

#Calculating the frequency could be done with angles but since it starts at  0, 0 let's add one round each time y starts to decrease again (See)
angles = np.angle((x - center[0]) + (y - center[1])*1.0j) 	#By using imaginary numbers we bypass all the stuff with quadrants and angles

# np.savetxt("foo.csv", angles, delimiter="\t")

#Now we need to check how many times the angles go back down to 0ish and add one to the rotationcounter
for i in range (1,steps -1):
	if np.abs(angles[i + 1] - angles[i]) > np.pi:
		rotations += 1
		lastRotationTime = i*h
		if startRotationTime == 0:	#starts clock
			startRotationTime = i*h
			rotations = 0

print rotations

# pl.figure()
# pl.plot(time, angles)		#This wasn't a really good plot for anything

pl.show()


print 'The center is located at '  + str(center)	#Just to make sure that enough gyrations is done for the mean method of finding the center is valid
print ' The theoretical gyration radius is ' + str(m*vx[0]/(q*B0))
print ' From the simulation it is ' + str(rho)

print ' The theoretical gyration frequency is ' + str(q*B0/(m))
print ' From the simulation it is ' + str(rotations / (lastRotationTime - startRotationTime) * 2.*np.pi)

