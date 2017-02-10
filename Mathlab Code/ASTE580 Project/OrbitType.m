function [ecc, orbit_type]=OrbitType(R, V, mu)
%% Description
% Name: OrbitPropagator
% By: Rob Holland
% Purpose: Determine an orbit's eccentricity and thus type given a 
% radius (R), velocity (V), and central body mu (mu).
%% Unit Notes
% e is unitless
% R0 is in km, V0 is in km/s
% Type is a string
%% Compute Eccentricity Vector sqrd (ECC2) and then Eccentricity (ecc)
r=sqrt(dot(R,R));
ECC2=1/mu*((dot(V,V)-mu/r)*R-dot(R,V)*V);
ecc=sqrt(dot(ECC2,ECC2));
%% Determines Orbit Type
if  ecc<1
    orbit_type=char('Elliptical');
elseif ecc>1
    orbit_type=char('Hyperbolic');
elseif ecc==0
    orbit_type=char('Circluar');
elseif ecc==1
    orbit_type=char('Parabolic');
end