function [a, i, raan, arg_p, theta]=ClassicElements(R, V, mu, ecc)
%% Description
% Name: ClassicElements
% By: Rob Holland
% Purpose: Computes the classical orbital elements given radius, velocity,
% eccentricity, and central body's mu.
%% Unit Notes
% p, a, r_p, & r_a are all in km
% i, raan, arg_p, theta are processed in radians and converted to
% degrees for display.
% e is unitless
% R is in km, V is in km/s
%% Unit Vectors
I=[1;0;0];
J=[0;1;0];
K=[0;0;1];
%% Computing a (km)
r=sqrt(dot(R,R));
v=sqrt(dot(V,V));
nrg=v^2/2-mu/r;
a=-mu/(2*nrg);
%% Computing i (rad)
h=cross(R,V);
h_hat=h/sqrt(dot(h,h));
i=acos(dot(h_hat,K));
%% Computing raan (rad)
N=cross(K,h_hat);
N_hat=N/sqrt(dot(N,N));
raan=acos(dot(N_hat,I));
%% Computing arg_p (rad)
ECC=1/mu*[(v^2-mu/r)*R-dot(R,V)*V];
arg_p=acos(dot(ECC,N_hat)/ecc);
if ECC(3) < 0
    arg_p=2*pi-arg_p;
end
%% Computing theta (rad)
theta=acos((a/r*(1-ecc^2)-1)/ecc);