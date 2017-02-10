%% Program Info
% Title: PatchedConic.m
% Written By: Rob Holland
% Purpose: Produce the following:
    % (1) required delta-V for injection
    % (2) Post-flyby values of r, v, & Beta
    % (3) Post-flyby values for helicentricellipse values of e, Theta, h, r_p,
    % a, & Tau

%% Data Clearing & Formating
clear
clc
format compact

%% Constants:
Astro_Unit=1.4958E8; % km, 1 AU for reference
    
%% Assumptions:
% All units are SI
% Departure Body is Earth/Terra
Mu_terra=398600.433; % km^3/s^2
r_terra=6378.14; % km (Earth Mean Equatoral Radius)
a_terra=1*Astro_Unit; %km (Semimajor Axis of Earth around the Sun)
% Primary Central Body is Sol
Mu_sol=132712440017.987; % km^3/s^2
    
%% Inputs:
% This section obtains the inputs from the user
% Parking Orbit Altitude in km from Departure Body Surface (h_parking)
h_parking=input('Enter the Initial Parking Orbit Altitude(km): ');
% Target Body Mu (mu_target)
Mu_target=input('Enter the Gm or Mu of the Target (km^3/s^2): ');
% Target Body Mean Radius (r_target)
r_target=input('Enter the Mean Surface Radius of the target (km): ');
% Target Body Semimajor Axis (a_target)
a_target=input('Enter the Semimajor Axis of the Target Heliocentric Orbit (AU): ');
a_target=a_target*Astro_Unit;
% Target Body Flyby Altitude in km
h_flyby=input('Enter the desired Flyby Altitude (km): ');
% Generalizing Assumptions for potential code re-use
a_initial=a_terra;  % Initial Body semimajor axis
r_initial=r_terra;  % Initial Body Mean Equatorial Radius
Mu_initial=Mu_terra; % Initial Body Mu
Mu_central=Mu_sol; % Mu of the major central body for the transit orbit
% Post input calculations
r_parking=h_parking+r_initial; % Parking orbit radius from Initial Body
r_flyby=h_flyby+r_target;   % Flyby radius from Target Body


%% Transit Orbit Values
% This section determines the properties for the transient orbit for use
% later. Hohmann transfer assumed for the transit orbit.
%Semimajor axis of the transit orbit (km)
a_transit=(a_initial+a_target)/2; 
% Circular Orbit Velocity of the Initial Body (km/s)
V_initial=sqrt(Mu_central/a_initial); 
% Circular Orbit Velocity of the Initial Body (km/s)
V_target=sqrt(Mu_central/a_target); 
% Transit Orbit Velocity @ Initial Body (km/s)
V_plus_initial=sqrt(Mu_central*(2/a_initial-1/a_transit)); 
% Transit Orbit Velocity @ Target Body (km/s)
V_minus=sqrt(Mu_central*(2/a_target-1/a_transit)); 

%% Injection Velocity
% This section determines the injection velocity required to achieve the 
% desired transit orbit (V_injection, km/s).
%Parking Orbit Velocity 
V_parking=sqrt(Mu_initial/r_parking);  % km/s
% Velocity at Infinity from Initial Body
V_inf_initial=V_plus_initial-V_initial;  % km/s
% Needed Velocity at the parking radius to acheive the Velocity at infinity
V_p_initial=sqrt(V_inf_initial^2+2*Mu_initial/r_parking); % km/s
% Injection delta V *OUTPUT*
V_injection=V_p_initial-V_parking % km/s

%% Hyperbolic Properties Post-Flyby at the Target Body
% Determines the various velocities and angles at the target body
% Infinite Velocity upon arrival at the target body (km/s)
V_inf_target=V_minus-V_target; 
% Hyperbolic eccentricity
e_hyper=1+r_flyby*V_inf_target^2/Mu_target;
% Simple Error Check!
if e_hyper<1
   error('Error! The Hyperbolic Eccentrcity cannot be less than 1')
end
% Delta hyperbolic angle (rad)
Delta=2*asin(1/e_hyper);
% Velocity departing the target body (km/s) *OUTPUT*
V_plus_target=sqrt(V_inf_target^2+V_target^2-2*V_inf_target*V_target*cos(pi-Delta));
% Determination of Beta departing the target body (rad) *OUTPUT*
Beta_target=acos((V_target^2+V_plus_target^2-V_inf_target^2)/(2*V_target*V_plus_target));
Betadeg=Beta_target*180/pi

%% Heliocentric Ellipse Properties
% Determines the eccentricy (e), Theta, momentum (h), radius perigee (r_p),
% semimajor axis (a), and period (Tau) for the ellipitcal heliocentric
% orbit following the flyby of the target.
% *OUTPUT* Data
R_central=a_target % km
V_central=V_plus_target % km/s
X0=R_central*V_central^2/Mu_central; 

%Eccentricty of new Heliocentric Orbit *OUTPUT*
e=sqrt((X0-1)^2*cos(Beta_target)^2+sin(Beta_target)^2)

%True Anomaly of new Heliocentric Orbit (radians) *OUTPUT*
Theta=atan((X0*sin(Beta_target)*cos(Beta_target))/(X0*cos(Beta_target)^2-1));
if Beta_target<0 || Beta_target>pi 
    Theta=Theta+pi;
else
    Theta;
end
Thetadeg=Theta*180/pi

%Momentum Magnitude of new Heliocentric Orbit (km^2/s^2) *OUTPUT*
h=R_central*V_central*cos(Beta_target) 

%Radius of Perigee of new Heliocentric Orbit (km) *OUTPUT*
r_p=h^2/Mu_central/(1+e*cos(Theta)) % km

%Semimajor Axis of new Heliocentric Orbit (km) *OUTPUT*
a=h^2/(Mu_central*(1-e^2))

%Period of new Heliocentric Orbit (s) *OUTPUT*
Tau=2*pi*sqrt(a^3/Mu_central)
