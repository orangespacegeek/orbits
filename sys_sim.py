import numpy as np
import datetime
import math

pi = math.pi

#Solar System Simulator
#Downloads the orbital parameters of the major planetary bodies and several large 
#asteroids and based on their orbital data, projects them forward to the indicated
#time.

#calculates the julian date
def julian(yr, mm, dd, ut):
	#calculates the julian date given the Georgian year (yr), month (mm), day (dd), and UTC (ut)
	JD = 367*yr-truncate(7*(yr+truncate((mm+9)/12)/4)) + truncate(275*mm/9) + dd + 1721013.5 + ut/24 - 0.5 * np.sign((100*yr+mm-190002.5)) + 0.5
	return JD


def truncate(n):
	# truncates the solution to the units digit
	m = n-n%1
	return m

def sys_julian():
	sys_t = datetime.datetime.utcnow()
	sys_hr = sys_t.hour+(sys_t.minute+sys_t.second/60)/60
	sys_julian=julian(sys_t.year, sys_t.month, sys_t.day,sys_hr)
	return sys_julian

def extract_TLE(filename, planet):
  """
  Given a file name for the planets (JPL p_elem_t2.txt), extract the two line
  element (TLE) set for the given planet
  Outputs TLE in the following format (a, e, i, L, lp, ln, da, de, di, dL, dlp, dln)
  """
  nxt_line = False
  TLE = []
  file = open(filename, 'rU')
  l_num = 0
  adjust = []

  for line in file:
  	l_num = l_num +1
  	if  line[:len(planet)] == planet and l_num>40 and planet == ('Jupiter' or 'Saturn' or 'Neptune' or 'Uranus' or 'Pluto'):
		#extracts adjustment for outer planets
  		adjust.append(float(line[10:22])) #extracts b
  		adjust.append(float(line[25:36])) #extracts c
  		adjust.append(math.radians(float(line[39:40]))) #extracts s
  		adjust.append(math.radians(float(line[53:64]))) #extracts f
  	elif line[:len(planet)] == planet and l_num<40:
  		TLE.append(float(line[10:21])) #extracts the semimajor axis (AU) from the TLE like those in the p_elem_t2.txt file
  		TLE.append(float(line[26:37])) #extracts the eccentricity (rad) from the TLE like those in p_elem_t2.txt file
  		TLE.append(math.radians(float(line[42:53]))) #extracts the inclination (deg) from the TLE like those in p_elem_t2.txt file, also converts to radians
  		TLE.append(math.radians(float(line[56:71]))) #extracts the L (deg) from the TLE like those in p_elem_t2.txt file, also converts to radians
  		TLE.append(math.radians(float(line[75:87]))) #extracts the longitude of perihelion (deg) from the TLE like those in p_elem_t2.txt file, also converts to radians
  		TLE.append(math.radians(float(line[90:103]))) #extracts the longitudal node (deg) from the TLE like those in p_elem_t2.txt file, also converts to radians
  		nxt_line = True
  	elif nxt_line:
  		TLE.append(float(line[10:21])) #extracts the semimajor axis rate (AU/Cy) from the TLE like those in the p_elem_t2.txt file
  		TLE.append(float(line[26:37])) #extracts the eccentricity rate (rad/Cy) from the TLE like those in p_elem_t2.txt file
  		TLE.append(math.radians(float(line[42:53]))) #extracts the inclination rate (deg/Cy) from the TLE like those in p_elem_t2.txt file, also converts to radians
  		TLE.append(math.radians(float(line[56:71]))) #extracts the L rate (deg/Cy) from the TLE like those in p_elem_t2.txt file, also converts to radians
  		TLE.append(math.radians(float(line[75:87]))) #extracts the longitude of perihelion rate (deg/Cy) from the TLE like those in p_elem_t2.txt file, also converts to radians
  		TLE.append(math.radians(float(line[91:103]))) #extracts the longitudal node rate (deg/Cy) from the TLE like those in p_elem_t2.txt file, also converts to radians
  		nxt_line = False

  	#Outputs TLE in the following format (a, e, i, L, lp, ln, da, de, di, dL, dlp, dln)
  	#Outputs the adjustment variables
  return (TLE, adjust)

def extract_planets():
  #Extracts all 8 planets of the Sol system from p_elem_t2.txt and puts them into a dictonary
  TLE_file = 'p_elem_t2.txt'
  planets = ('Mercury', 'Venus', 'EM Bary', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune')
  TLE_DB = {}
  for planet in planets:
  	TLE_DB[planet] = extract_TLE(TLE_file, planet)
  return TLE_DB


def adv_planet(TLE, date, adjust):
  #advances individual body from p_elem_t2.txt and outputs an updated TLE
  T = (date-2451545)/36525 #Time adjustment for equations to use best fit
  # updates the TLE to the appropriate position based on the Julian date
  a = TLE[0]+TLE[6]*T
  e = TLE[1]+TLE[7]*T
  i = TLE[2]+TLE[8]*T
  L = TLE[3]+TLE[9]*T
  longperi = TLE[4]+TLE[10]*T
  longnode = TLE[5]+TLE[11]*T
  # checks the adjust input and assigns b,c,s,f if available or zeros 
  if len(adjust) > 0:
  	b = adjust[0]
  	c = adjust[1]
  	s = ajdust[2]
  	f = adjust[3]
  else:
  	b = 0
  	c = 0
  	s = 0
  	f = 0
  # computes the additional traditional TLE for R & V computation, 
  # argument of perihelion (arg_peri) and the mean anomaly (mean)
  arg_peri = longperi - longnode
  mean_anom = L - longperi + b*T*T+c*np.cos(f*T) + np.sin(f*T)
  
  M_0 = mean_anom%pi
  E = M_0 + e*np.sin(M_0)
  M = M_0
  dE = 1
  n = 0
  while dE>0.0000001 and n < 1000000:
  	dM = M - (E - e*np.sin(E))
  	dE = dM/(1-e*np.cos(E))
  	E = E + dE
  	M = E - e*np.sin(E)
  	n = n + 1

  TLE_nu = [a, e, i, longnode, arg_peri, E]

  return TLE_nu

def radius(TLE):
	#Computes the radius vector of the object given its current TLE to a reference plane
	#For heliocentric objects, reference plane is J2000 ecliptic
	#Heliocentric coordinates in its orbital plane, x-axis aligned from focus to perihelion
	a = TLE[0]
	ecc = TLE[1]
	I = TLE[2]
	longnode = TLE[3]
	arg_peri = TLE[4]
	E = TLE[5]
	x_prime = a*(np.cos(E)-ecc)
	y_prime = a*np.sqrt(1-ecc*ecc)*np.sin(E)
	z_prime = 0

	x_ecl = float((np.cos(arg_peri)*np.cos(longnode)-np.sin(arg_peri)*np.sin(longnode)*np.cos(I))*x_prime + (-np.sin(arg_peri)*np.cos(longnode)-np.cos(arg_peri)*np.sin(longnode)*np.cos(I))*y_prime)
	y_ecl = float((np.cos(arg_peri)*np.sin(longnode)+np.sin(arg_peri)*np.cos(longnode)*np.cos(I))*x_prime + (-np.sin(arg_peri)*np.sin(longnode)+np.cos(arg_peri)*np.cos(longnode)*np.cos(I))*y_prime)
	z_ecl = float((np.sin(arg_peri)*np.sin(I))*x_prime + (np.cos(arg_peri)*np.sin(I))*y_prime)

	R_vec = [x_ecl, y_ecl, z_ecl]

	return R_vec

def velocity(TLE):
	#Computes the velocity vector of the object given its current TLE
	V_vec = 0
	return V_vec

def main():
	
	JD_now = sys_julian()
	
	TLE_DB = extract_planets()

	TLE_EM_now = adv_planet(TLE_DB['EM Bary'][0], JD_now, TLE_DB['EM Bary'][1])

	TLE_Mars_now = adv_planet(TLE_DB['Mars'][0], JD_now, TLE_DB['Mars'][1])

	R_EM = radius(TLE_EM_now)

	R_Mars = radius(TLE_Mars_now)

	dist = vec_sep(R_Mars, R_EM)

	print ('')

	print('Earth-Moon TLE', TLE_EM_now)
	print('Earth-Moon Semi-Major Axis (AU):', TLE_EM_now[0]) #prints semi-major axis for Earth-Moon
	print('Earth-Moon Radius (AU): ', np.sqrt(R_EM[0]*R_EM[0]+R_EM[1]*R_EM[1]+R_EM[2]*R_EM[2])) #prints updated Earth-Moon heliocentric radius
	
	print('Mars Semi-Major Axis (AU):', TLE_Mars_now[0]) #prints semi-major axis for Earth-Moon
	print('Mars Radius (AU):', np.sqrt(R_Mars[0]*R_Mars[0]+R_Mars[1]*R_Mars[1]+R_Mars[2]*R_Mars[2])) #prints updated Mars heliocentric radius
  #spacer
  
	print('Seperation (AU):', dist)
	print ('')

if __name__ == '__main__':
  main()
