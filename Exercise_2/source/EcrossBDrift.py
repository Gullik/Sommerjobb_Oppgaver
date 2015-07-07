import numpy as np
import scipy as sp
import pylab as pl
import scipy.stats as stats



#Declaration of several variables, all the units are in SI
e =  -1.602E-19		#C
m_e = 9.109E-31		#kg
m_i = 2*1.66E-27 	#kg
B0 =  50000E-9		#Tesla


particle = "electron"
# particle = "ion"

# case = "ExB"
case = "VaryV"

E_x = 50 * 10**-3		#mV/m
E_y = 0

# #Choose if for ion or electron
if particle == "electron":
	# #electron settings
	h = 1.E-14		#stepsize
	m = m_e	
	q = e

if particle == "ion":
	#Ion settings
	h = 1.E-8
	m = m_i
	q = 2*e

gyrations = 4

steps = int( np.abs(m/(q*B0)/h) * 2*np.pi  * gyrations)	#Since the gyration has a frequency qB/m, it's period is m/(qB) 


#Defining the vectors which contains the positions and velocities during the simulation
x = np.zeros(steps)
y = np.zeros(steps)
vx = np.zeros(steps)
vy = np.zeros(steps)
kineticEnergy = np.zeros(steps)
angles = np.zeros(steps)
time = np.linspace(0,steps*h, steps) #Only used for plotting the kinetic energy

#Rewriting it to use vectors instead of x, y, vx and vy
# r = n




# pl.figure()
# pl.plot(time, kineticEnergy)		#This wasn't a really good plot for anything

# pl.show()

if case == "ExB":
	#Setting initial conditions
	X = 0
	Y = 0
	# VX = 500 + E_y/B0 
	# VY = -E_x /B0
	velocityMagnitude = 500

	angle = np.pi

	VX = np.sin(angle) * velocityMagnitude
	VY = np.cos(angle) * velocityMagnitude 


	x[0] = X
	y[0] = Y
	vx[0] = VX
	vy[0] = VY

	#Precalculating some constants so not to calculate them in the loop
	C1 = h*q/m
	C2 = h*(q*B0/m)

	for i in range(0,steps -1):
		#Rewrote it so it doesn't go searching through the vectors so often, faster

		#Let us first calculate the steps since last move
		dX = h * VX
		dY = h * VY
		dVX = C1*E_x + C2*VY
		dVY = C1*E_y - C2*VX

		#Calculate new positions and velocities
		X = X + dX
		Y = Y + dY
		VX = VX + dVX
		VY = VY + dVY

		#Storing datapoints
		vx[i + 1] = VX
		vy[i + 1] = VY
		x[i + 1] = X
		y[i + 1] = Y
		if i % int(0.25*steps) == 0:
			print str( float(i) / float(steps) * 100)  + '%'


	#Calculating the kinetic energy to see that it stays decently constant
	kineticEnergy = m*(vx*vx + vy*vy)/2.
	print "relative change E_K = " + str((np.max(kineticEnergy) - np.min(kineticEnergy))/np.mean(kineticEnergy)) 
	print "std(E_K) = " + str(stats.tstd(kineticEnergy))
	pl.figure()
	pl.plot (time, kineticEnergy)
	pl.xlabel('s')
	pl.ylabel('J')
	pl.title('Kinetic Energy')
	pl.savefig('kineticEnergy.eps')

	# pl.figure()
	f, axArray = pl.subplots(1,2 ,sharex = True, sharey = True)
	pl.suptitle('Gyration of an ' + str(particle))

	axArray[0].plot(x,y)
	axArray[0].set_title('Numerical')
	axArray[0].set_xlabel("x [m]")
	axArray[0].set_ylabel("y [m]")


	#Calculate the theoretical drift
	v_D = -E_x/B0	#This is in the y direction, it is zero in all the other directions

	time = np.linspace(0,steps*h,1000) #Only used for plotting the kinetic energy
	omega_c = (q*B0/m)
	rho_c = np.abs(m*np.sqrt( (vx[0]+ E_y/B0 )*(vx[0] + E_y/B0 ) + (vy[0] - E_y/B0 )*(vy[0] - E_y/B0 ) )/(q*B0)) # vx[0]*vx[0] + vy[0]*vy[0] )/( q*B0 ))
	phase = np.pi/2.

	print np.sqrt( (vx[0] + E_y/B0)*(vx[0] + E_y/B0) + (vy[0] - E_x/B0)*(vy[0] - E_x/B0)) 

	xAnalytical = np.sign(q)*(rho_c*np.sin(phase) -  rho_c*np.sin( omega_c*time + phase )) 
	yAnalytical = np.sign(q)*(rho_c*np.cos(phase) -  rho_c*np.cos( omega_c*time + phase )) + v_D*time


	axArray[1].plot(xAnalytical,yAnalytical)
	# axArray[1].text(xAnalytical[0], yAnalytical[0], 'start', fontsize=15)
	axArray[1].set_title('Analytical')
	axArray[1].set_xlabel("x [m]")
	f.savefig("ExB" + particle + ".eps")
	pl.show()

	#Calculating the numerical ExB drift, I will just do it after an integer number of gyrations so the position of the particle in the gyration is not important
	tEnd = np.abs(m/(q*B0)) * 2*np.pi  * gyrations	#The period times number of wanted gyrations
	tStart = 0

	driftx = x[steps-1] - x[0] / tEnd	#This should be the movement after a certain number of gyrations
	drifty = y[steps-1] - y[0] / tEnd
	
	print 'Numerical Drift in x direction: ' + str(driftx)
	print 'Analytical drift in x direction: ' + str(E_y/B0)

	print 'Numerical Drift in y direction: ' + str(drifty)
	print 'Analytical drift in y direction: ' + str(-E_x/B0*tEnd)


if case == "VaryV":

	#Setting up subplots
	f, axArray = pl.subplots(1,4 ,sharex = True, sharey = True)
	pl.suptitle('Gyration of an ' + str(particle))

	axArray[0].set_xlabel("x [m]")

	for i in range(0,4):
		#Setting initial conditions
		X = 0
		Y = 0
		# VX = 500 + E_y/B0 
		# VY = -E_x /B0
		velocityMagnitude = 500

		angle = i*np.pi/2.

		VX = np.sin(angle) * velocityMagnitude
		VY = np.cos(angle) * velocityMagnitude 


		x[0] = X
		y[0] = Y
		vx[0] = VX
		vy[0] = VY

		#Precalculating some constants so not to calculate them in the loop
		C1 = h*q/m
		C2 = h*(q*B0/m)

		for i in range(0,steps -1):
			#Rewrote it so it doesn't go searching through the vectors so often, faster

			#Let us first calculate the steps since last move
			dX = h * VX
			dY = h * VY
			dVX = C1*E_x + C2*VY
			dVY = C1*E_y - C2*VX

			#Calculate new positions and velocities
			X = X + dX
			Y = Y + dY
			VX = VX + dVX
			VY = VY + dVY

			#Storing datapoints
			vx[i + 1] = VX
			vy[i + 1] = VY
			x[i + 1] = X
			y[i + 1] = Y
			if i % int(0.25*steps) == 0:
				print str( float(i) / float(steps) * 100)  + '%'

		# pl.figure()

		axArray[i].set_ylabel("y [m]")


		#Calculate the theoretical drift
		v_D = -E_x/B0	#This is in the y direction, it is zero in all the other directions

		time = np.linspace(0,steps*h,1000) #Only used for plotting the kinetic energy
		omega_c = (q*B0/m)
		rho_c = np.abs(m*np.sqrt( (vx[0]+ E_y/B0 )*(vx[0] + E_y/B0 ) + (vy[0] - E_y/B0 )*(vy[0] - E_y/B0 ) )/(q*B0)) # vx[0]*vx[0] + vy[0]*vy[0] )/( q*B0 ))
		phase = np.pi/2.

		print np.sqrt( (vx[0] + E_y/B0)*(vx[0] + E_y/B0) + (vy[0] - E_x/B0)*(vy[0] - E_x/B0)) 

		xAnalytical = np.sign(q)*(rho_c*np.sin(phase) -  rho_c*np.sin( omega_c*time + phase )) 
		yAnalytical = np.sign(q)*(rho_c*np.cos(phase) -  rho_c*np.cos( omega_c*time + phase )) + v_D*time



		f.savefig("ExBVaryingV" + particle + ".eps")
		pl.show()

	