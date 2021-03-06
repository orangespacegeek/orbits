function dm=PropMassReq(dV,m1,Isp)
% Computes the PropMass expelled given dV required (km/s), mass initial 
%(kg), and system specific impulse (s).
g0=9.81; %acceleration due to gravity (Earth) m/s^2
dV=dV*1000; %km/s to m/s
% computes mass after the burn
m0=m1*exp(dV/(g0*Isp));
% determines mass burned in kg
dm=m0-m1;