import numpy as np
import pylab as plt
import sys as sys

def epsilon(omega, k):
	cs	= 1000 				#m/s 		Speed of sound
	O_i = 300				#1/s 		Ion gyrating frequency
	O_e = 8.E6 				#1/s 		Electron gyrating frequency
	theta = np.deg2rad(10)	#rad 		angle between propagation and B
	tan2theta = np.tan(theta)**2

	# print 'Second term: ' + str(np.min((k*k*cs*cs)/(omega*omega)))

	# print 'Third term: ' + str (O_i*tan2theta/omega*(
	# 				1/(	omega*tan2theta/O_e  -  O_e/omega *(1 - omega*omega/O_e) ) -
	# 				1/( omega*tan2theta/O_i  -  O_i/omega *(1 - omega*omega/O_i) )  ) )

	# print str(O_i*tan2theta/omega*(
	# 				1/(	omega*tan2theta/O_e  -  O_e/omega *(1 - omega*omega/O_e) ) -
	# 				1/( omega*tan2theta/O_i  -  O_i/omega *(1 - omega*omega/O_i) )  ) /
	# 				np.min((k*k*cs*cs)/(omega*omega)))

	# print 'Second term: ' + str( 1/(	omega*tan2theta/O_e  -  O_e/omega *(1 - omega*omega/O_e) ) )

	# print 'Third term: ' + str( 1/( omega*tan2theta/O_i  -  O_i/omega *(1 - omega*omega/O_i) ) )


	return 1 - (k*k*cs*cs)/(omega*omega) + O_i*tan2theta/omega*(
					1/(	omega*tan2theta/O_e  -  O_e/omega *(1 - omega*omega/O_e) ) -
					1/( omega*tan2theta/O_i  -  O_i/omega *(1 - omega*omega/O_i) )  )


def epsilon2(omega, k):
	lam_se = 1
	O_e = 1

	return 1 + 5./3.*k*k - omega*omega/(O_e*O_e)


def bisection(function, k, xmin, xmax, iterations):
	
	fmin = function(xmin , k)
	fmax = function(xmax , k)

	xx = np.arange(xmin, xmax, (xmax-xmin)/1000)

	if fmin*fmax > 0:
		print '0 or more than 1 bisection in the interval'
		return 0

	if fmin > 0:
		#If it is decreasing swap all the values
		temp = xmin
		xmin = xmax
		xmax = temp

	# print "Got here"
	dx = (xmax - xmin)*0.5
	for i in range(iterations):
		
		xmid = xmin + dx

		if function(xmid, k) < 0:
			xmin = xmid
		else:
			xmax = xmid
		dx *= 0.5

		# plt.figure()
		# plt.plot(xx, function(xx, k), label = '$\epsilon(\omega,k)$')
		# plt.grid()
		# plt.axvspan(xmin, xmax, facecolor='g', alpha=0.2, label = 'root interval')
		# plt.axhline( y = 0, color = 'black')
		# plt.axvline( x = (xmax + xmin)*0.5, color = 'black' )
		# plt.legend()
		# plt.xlabel('$\omega$')
		# plt.ylabel('$\epsilon$')
		# plt.savefig('../figures/bisection_' + str(k) + '_' + str(i) + '.eps')
		# plt.show()



	plt.close('all')
	
	return (xmax + xmin)*0.5, np.abs(dx)


if __name__ == '__main__':

	# xmin =  np.power(10,int(sys.argv[1]))
	# xmax =  np.power(10,int(sys.argv[2]))
	# xmin = -1.E3
	xmin =  0.00001
	xmax =  np.power(10,6)

	kk = [1, 10, 200, 300, 400, 500, 600, 700, 800, 900, 999]

	# kk = np.arange(0,3, 0.1)
	# kk = [100]
	# kk = [1, 2, 3 , 4, 5, 10, 20, 30 ,50]#, 10000, 1000000]
	# kk[0] = 1
	roots = np.zeros(len(kk))

	i = 0
	for k in kk:
		iterations = 50
		print 'Finding root with k = ' + str(k)

		omega = np.arange(0.0001,1.E6, 10)

		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.plot(omega, epsilon(omega,k))
		ax.set_ylim(-10, 10)
		ax.set_xscale('log')
		plt.axhline( y = 0, color = 'black')
		plt.grid()

		plt.show()

		xmin = float( pow(10, int(raw_input("xmin: "))) )
		xmax = float( pow(10, int(raw_input("xmax: "))) )

		print 'Root is estimated to be in the interval: (' + str(xmin) + ',' + str(xmax) + ')'

		# roots[i], domega = bisection(epsilon, k, xmin, xmax, iterations)

		print roots[i]
		i += 1

		# plt.show()

	plt.close('all')
	plt.figure()
	plt.plot(kk, roots)
	plt.xlim([0, kk[-1]])
	plt.xlabel('$k/\lambda_{se}$')
	plt.ylabel('$ \omega/\omega_{pe} $')
	plt.title('Dispersion diagram for electron plasma waves')
	plt.ylim([0, roots[-1]])
	# plt.savefig('../figures/simple_dispersion.eps')
	plt.savefig('../figures/ionAcousticWaves')
	plt.show()

