% ASTE580 20113     Holland, Robert
% HW#5 3-1d
% Rework 3-1a, b, & c in Matlab
clear
clc
format compact

% Constants
mu_e=398600; % km^3/s^2

% Initial Conditions Query
v_0=9;       % km/s
fprintf('The initial velocity (v_0) is %1.1f km/s.\n', v_0); 
r_0=7500;    % km
fprintf('The initial radius (r_0) is %4.1f km.\n', r_0); 
B_0=0;       % deg/radians
fprintf('The initial flight path angle (B_0) is %4.1f deg.\n', B_0); 

% a) Circular orbit?
disp('a) Is this orbit circular?')

% To test for a circular orbit, determine whether the energy (nrg_0) at t0 is equal
% to the energy of a circular orbit (nrg_C) of the same radius.  The radius
% (r_C) and semimajor axis (a_C) of a circle are the same.

nrg_0=v_0^2/2-mu_e/r_0;
a_C=r_0;
nrg_C=-mu_e/(2*a_C);
fprintf('- The energy of this orbit is %2.2f km^2/s^2.\n', nrg_0);
fprintf('- The energy of a circular orbit, with a radius of %4.1f km, is %2.2f km^2/s^2.\n', a_C, nrg_C); 
if nrg_0==nrg_C disp('Ans: Since energy is equal, the the orbit is circular.');
    else disp('Ans: The two energies are not equal, therefore the orbit is not circular.');
end

disp('b) Where is the spacecraft in its orbit?');
fprintf('-A flight path angle of %3.2f indicates the spacecraft is either at apogee or perigee\n', B_0);
% To test for the vehicle's location we need to determine the semimajor
% axis of the spacecraft and whether it is greater than or less than the
% given radius.

a_0=-mu_e/(2*nrg_0);
fprintf('-The semimajor axis is %3.2f km\n', a_0);
if r_0<a_0 disp('Ans: The spacecraft is at perigee');
    r_p=r_0;
    e=1-r_p/a_0;
elseif r_0>a_0 disp('Ans: The spacecraft is at apogee');
    r_a=r_0;
    e=r_a/a_0-1;
end

disp('c) What are the values for energy (km^2/s^2), a (km), p (km), and h (s^2/km^2)')

energy=nrg_0
a=a_0
p=a*(1-e^2)
h=sqrt(p*mu_e)

