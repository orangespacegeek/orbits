import sys
import numpy as np
import datetime
import math

pi = math.pi

#Keplar Functions
#Provides functions for orbital mechanics based on Keplar's laws

def truncate(n):
	# truncates the solution to the units digit
	m = n-n%1
	return m

def aperiod(mu,period):
	# solves for the semimajor axis 
	 return (mu*(period/(2*pi))**2)**(1/3)

def ClassicElements(R,V,mu):
	# Computes the classical orbital elements given radius, velocity, and
	# central body's mu.
	# Unit Notes
	# p, a, r_p, & r_a are all in km
	# i, raan, arg_p, theta are processed in radians and converted to degrees for display.
	# ecc is unitless
	# R is in km, V is in km/s
	# Unit Vectors
	I=np.array([1,0,0])
	J=np.array([0,1,0])
	K=np.array([0,0,1])
	
	# Computing a (km)
	r=np.sqrt(np.dot(R,R))
	v=np.sqrt(np.dot(V,V))

	nrg=(v**2)/2-mu/r
	print(nrg)
	a=-mu/(2*nrg)
	print(a)
	# Computing ecc (unitless)
	H = np.cross(R,V)
	# ECC=((v**2-mu/r)*R-np.dot(R,V)*V)/mu
	ECC = np.cross(V,H)/mu-R/r
	ecc=np.linalg.norm(ECC)

	# Computing i (rad)
	h_hat=H[2]/np.sqrt(np.dot(H,H))
	i=np.arccos(h_hat)
	if i > pi:
		i = 2*pi-i

	# Computing raan (rad)
	N=np.cross(K,H)
	N_mag = np.sqrt(np.dot(N,N))
	N_hat=N[0]/N_mag
	raan=np.arccos(N_hat)
	if N[1] < 0: #WARN these seem to be opposite of the normal approach
		raan = 2*pi-raan

	# Computing arg_p (rad)
	arg_p=np.arccos(np.dot(ECC,N)/(ecc*N_mag))
	if ECC[2] < 0: #WARN these seem to be opposite of the normal approach
		arg_p=2*pi-arg_p

	# Computing theta (rad)
	theta=np.arccos((a/r*(1-ecc**2)-1)/ecc)
	if np.dot(R,V) < 0 and theta > pi:
		theta = 2*pi - theta

	# Packages the elements for output
	TLE_classic = [a, ecc, i, raan, arg_p, theta]
	
	return TLE_classic

def propagate(TLE, mu, dt):
	# Given an initial radius (R0) & initial velocity (V0) for an orbit, central body, 
	# and the amount of time to propagate to a new radius (Rn) & new velcoity (Vn).
	# Unit Notes
	# mu is in km^3/s^2
	# R0, R_cb & Rn are in km, V0, Vn, V_inf are in km/s
	# p, a, r_p, & r_a are all in km
	# i, raan, arg_p, theta, E, M are processed in radians
	# dt, period, and t_p (time to perigee) are processed in seconds
	# ecc is unitless.
	a = TLE[0]
	ecc = TLE[1]
	i = TLE[2]
	raan = TLE[3]
	arg_peri = TLE[4]
	theta = TLE[5]
	# Determine the period and the time remaining after several orbits
	
	# Check if orbit is elliptical or hyperpolic and determines E or F
	if ecc < 1: #elliptical
		E = np.arccos((ecc + np.cos(theta))/(1+ecc*np.cos(theta))) 
		TP = 2*pi*np.sqrt(a**3/mu)
		dt_p = dt%TP
		if dt_p == 0:
			return TLE
		else:
			a = a
			#continue with program

		#Netwons Method
		M1  = np.sqrt(mu/a**3)*dt_p
		while M1 - E > 10**-6:
			M2 = E - ecc*np.sin(E)
			E = E + (M1-M2)/(1 - ecc*np.cos(E))
			M1 = M2
		#Solve for the new true anomaly (theta_nu)
		top = ecc - np.cos(E)
		bottom = ecc*np.cos(E)-1
		theta_nu = np.arccos(top/bottom)
	else:  #parabolic/hyperbolic
		F = np.arccosh((ecc + np.cos(theta))/(1+ecc*np.cos(theta)))
		print('WARN: Eccentricty is greater than  or equal to 1, orbit is parabolic or hyperbolic')
		return


	TLE_out = [a, ecc, i, raan, arg_peri, theta_nu]
	return TLE_out

def calc_vectors(TLE, mu):
	#Calculates the radius and velocity vectors given a TLE (a, ecc, i, raan, arg_p, theta)
	# mu is in km^3/s^2
	# p, a, r_p, & r_a are all in km
	# i, raan, arg_p, theta, E, M are processed in radians
	# dt, period, and t_p (time to perigee) are processed in seconds
	# ecc is unitless.
	# R in km, V in km/s
	a = TLE[0]
	ecc = TLE[1]
	i = TLE[2]
	raan = TLE[3]
	arg_peri = TLE[4]
	theta = TLE[5]

	#Solve for P
	p = a*(1-ecc**2)

	#Solve for the radius based on the TLE in the perifocal system
	r = p / (1+ecc*np.cos(theta))
	R_P = r*np.cos(theta)
	R_Q = r*np.sin(theta)
	R_W = 0
	R_peri = np.array([R_P,R_Q, R_W])
	
	#Solve for the velocity based on the TLE in the perifocal system
	V_P = np.sqrt(mu/p)*(-np.sin(theta))
	V_Q = np.sqrt(mu/p)*(ecc+np.cos(theta))
	V_W = 0
	V_peri = np.array([V_P,V_Q, V_W])
	
	#Rotational Matrix for Perifocal Coordinates (PQW) to interial ref frame (IJK)

	R11 = np.cos(raan)*np.cos(arg_peri)-np.sin(raan)*np.sin(arg_peri)*np.cos(i)
	R12 = -np.cos(raan)*np.sin(arg_peri)-np.sin(raan)*np.cos(arg_peri)*np.cos(i)
	R13 = np.sin(raan)*np.sin(i)

	R21 = np.sin(raan)*np.cos(arg_peri)+np.cos(raan)*np.sin(arg_peri)*np.cos(i)
	R22 = -np.sin(raan)*np.sin(arg_peri)+np.cos(raan)*np.cos(arg_peri)*np.cos(i)
	R23 = -np.cos(raan)*np.sin(i)

	R31 = np.sin(arg_peri)*np.sin(i)
	R32 = np.cos(arg_peri)*np.sin(i)
	R33 = np.cos(i)

	Rot = np.array([[R11,R12,R13],[R21,R22,R23],[R31,R32,R33]])

	#Rotate R_peri and V_peri
	R = np.matmul(Rot,R_peri)
	V = np.matmul(Rot,V_peri)

	#Output the R & V vectors
	output = [R, V]

	return output

def test(case):
	#Test the functionality of keplar
	#Case 1
	print('Case #:', case)
	if case == '1':
		mu=398600.433 #Earth
		R_obj=6378.14 #Earth
		R0=np.array([-14192.498, -16471.197, -1611.2886])
		V0=np.array([-4.0072937, -1.2757932, 1.9314620])
		dt=8*60*60 #8 hours

	elif case == '2':
		mu=132712440017.987
		R_obj=6.96E5
		R0=np.array([148204590.0357, 250341849.5862, 72221948.8400])
		V0=np.array([-20.5065125006, 7.8793469985, 20.0718337416])
		dt=((10*24+0)*60+0)*60+0
	elif case == '3':
		mu=37940626.061
		R_obj=60268
		R0=np.array([-321601.0957, -584995.9962,-78062.5449])
		V0=np.array([8.57101142, 7.92783797, 1.90640217])
		dt=((0*24+10)*60+47)*60+39.30;
	else:
		return		

	TLE1_0 = ClassicElements(R0,V0,mu)
	[R_chk, V_chk] = calc_vectors(TLE1_0,mu)
	print('R0:', R0)
	print('R_chk:', R_chk)
	print('V0:', V0)
	print('V_chk:', V_chk)
	TLE1_8 = propagate(TLE1_0,mu,dt)
	[R_8, V_8] = calc_vectors(TLE1_8, mu)
	print(TLE1_0[5]*180/pi)
	print(TLE1_8[5]*180/pi)

	return

def main():
	args = sys.argv[1:]

	if not args:
		print('usage: [test] [case] to test that keplar is functioning, otherwise, it is a funcction library for other programs')
		sys.exit(1)
	if args[0] == 'test':
		if len(args)<2:
			print('No Case given')
		else:
			case = args[1]
			test(case)
	else:
		print('usage: [test] [case] to test that sys_sim is functioning, otherwise, it is a funcction library for other programs')
  
if __name__ == '__main__':
  main()
