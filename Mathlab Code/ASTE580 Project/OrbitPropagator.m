%% Description
% Name: OrbitPropagator
% By: Rob Holland
% Purpose: Given an initial radius (R0) & initial velocity (V0) for an 
% orbit, central body, and the amount of time to propagate to a new 
% radius (Rn) & new velcoity (Vn).
%% Unit Notes
% mu is in km^3/s^2
% R0, R_cb & Rn are in km, V0, Vn, V_inf are in km/s
% p, a, r_p, & r_a are all in km
% i, raan, arg_p, theta, E, M are processed in radians and converted to
% degrees for display.
% dt, period, and t_p (time to perigee) are processed in seconds and displayed
% in days, hrs, mins, and sec (dd hh:mm:ss).
% e is unitless.
%% Input & Orbit Type
% [R0, V0, mu, R_cb, dt]=UserInput; % Prompts user for input
tol=1E-6; % Tolerance value
[R0, V0, mu, R_cb, dt]=CaseInput; % Prompts user for case # to provide inputs
[ecc, orbit_type]=OrbitType(R0, V0, mu); %computes eccentricity given radius, velocity, and mu
if (ecc>1-tol) && (ecc<1+tol)
    [p, r_p, i, raan, arg_p, theta]=ParabolicElements(R0, V0, mu, ecc);
elseif ecc>1
    [p, a, V_inf, r_p, i, raan, arg_p, theta]=HyperbolicElements(R0, V0, mu, ecc);
elseif ecc<1
    [a, i, raan, arg_p, theta]=ClassicElements(R0, V0, mu, ecc);
end