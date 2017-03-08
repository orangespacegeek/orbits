import numpy as np
import datetime

#---Solar System Constants---
#dictonary of Solar System gravity parameters (m^3/s^2)
#Lists taken from conslidation at: https://en.wikipedia.org/wiki/Standard_gravitational_parameter
#
mu_dict = {}
mu_dict['Sun'] = 1.327124400189*10**20
mu_dict['Mercury'] = 2.20329*10**13
mu_dict['Venus'] = 3.248599*10**14
mu_dict['Earth'] = 3.9860044189*10**14
mu_dict['Terra'] = mu_dict['Earth']
mu_dict['Mars'] = 4.2828372*10**13
mu_dict['Ceres'] = 6.26325*10**10
mu_dict['Jupiter'] = 1.266865349*10**17
mu_dict['Saturn'] = 3.79311879*10**16
mu_dict['Uranus'] = 5.7939399*10**15
mu_dict['Neptune'] = 6.8365299*10**15
mu_dict['Pluto'] = 8.719*10**11
mu_dict['Eris'] = 1.1089*10**12
mu_dict['Luna'] = 4.90486950*10**12
mu_dict['Moon'] = mu_dict['Luna']

#---Object/Vehicle Data---
# Initial Conditions for object of interest
body_0 = "Sun" #primary body that is being orbited at object
r_0 = np.array([1,2,3]) #Object radius vector relative to body_0 at t=0 (m)
v_0 = np.array([1,2,3]) #Object velocity vector relative to body_0 at t=0 (m/s)
#Geogorian time in Earth UTC
year_AD = 2017 #Georgorian year, the murder of Caesar August occuring in -44
month_Terra = 6 #numerical month, Jan = 1, Feb = 2, ... Nov = 11, Dec =12
day_Terra = 20 #numerical day of the month
hour_UT = 12 #hours in Universal Time, ie Zulu Time
min_UT = 0 #minutes in Universal Time, ie Zulu Time
sec_UT = 0 #seconds in Universal Time, ie Zulu Time
year_HE = year_AD + 10000 #conversion to the Human Era year for display purposes
UT = hour_UT+(min_UT+sec_UT/60)/60



#Assigns the mu for the body at t=0
mu_0 = mu_dict[body_0]

#calculates the julian date
def julian(yr, mm, dd, ut):
	JD = 367*yr-truncate(7*(yr+truncate((mm+9)/12)/4)) + truncate(275*mm/9) + dd + 1721013.5 + ut/24 - 0.5 * np.sign((100*yr+mm-190002.5)) + 0.5
	return JD

def truncate(n):
	m = n-n%1
	return m

def intial_dates_display(year, month, day, time):
	sys_t = datetime.datetime.utcnow()
	sys_hr = sys_t.hour+(sys_t.minute+sys_t.second/60)/60
	sys_julian=julian(sys_t.year, sys_t.month, sys_t.day,sys_hr)
	in_julian=julian(year, month, day, time)
	print ('System Julian Time:', sys_julian)
	print ('Trajectory Initialization Time:', in_julian)

	return in_julian

def main():
	
	print ('')
	print ('---Holland Trajectory Calculator---')
	t_julian0 = intial_dates_display(year_AD,month_Terra,day_Terra,UT)

	print ('')
if __name__ == '__main__':
  main()