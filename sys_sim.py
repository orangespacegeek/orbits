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
  a = TLE[0]
  e = TLE[1]
  i = TLE[2]
  L = TLE[3]
  longperi = TLE[4]
  longnode = TLE[5]
  da = TLE[6]
  de = TLE[7]
  di = TLE[8]
  dL = TLE[9]
  dlongperi = TLE[10]
  dlongnode = TLE[11]
  
  #for element in rankname_list:
  #	element_nu = element[4:-5]
  #	rank_name_nu.append(element_nu)
  #return rank_name_nu
  return

def main():
	
	JD_now = sys_julian()
	
	TLE_DB = extract_planets()

	TLE_now = adv_planet(TLE_DB['Mars'][0], JD_now, TLE_DB['Mars'][1])

	print ('')

	print('Mars: ', TLE_DB['Mars'][0]) #prints Mars adjustment []
	print('Jupiter:', TLE_DB['Jupiter'][0]) #prints Jupiter TLE

	print ('')

if __name__ == '__main__':
  main()
