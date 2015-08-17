import numpy as np
import pylab as plt

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
