function [R0, V0, mu, R_obj, dt] = UserInput
%% Description
% Name: UserInput
% By: Rob Holland
% Purpose: Ask the user to input the R0 & V0 vectors, central body, and
% propagation time for use elsewhere in the utilizing program.
% R0 & R_obj are in km, V0 is in km/s, mu is in km^3/s^2, dt is in s
%% User Request
fprintf('Central Body can be any of the Sun, the 8 Planets, the Moon, or Titan. \n');
CenBody=input('What is the central body for this Orbit?  ', 's');
fprintf('For R0 & V0 provide them in the following format [i;j;k]. \n');
fprintf('Please provide distance in km & velocity in km/s. \n');
R0=input('Provide the initial radius (R0) for the orbit:  ');
V0=input('Provide the initial velocity (V0) for the orbit:  ');
fprintf('Please provide the propagation time in seconds \n');
dt=input('Provide the duration of time (dt) to propagate:  ');
%% Setting mu & R_obj
% The following if loop assigns mu & R_obj based on which central body they
% input.
if strcmp(CenBody, 'Sun')==1
    mu=132712440017.987;
    R_obj=6.96E5;
elseif strcmp(CenBody, 'Mercury')==1
    mu=22032.08;
    R_obj=2439.7;
elseif strcmp(CenBody, 'Venus')==1
    mu=324858.599;
    R_obj=6051.8;
elseif strcmp(CenBody, 'Earth')==1
    mu=398600.433;
    R_obj=6378.14;
elseif strcmp(CenBody, 'Moon')==1
    mu=4902.801;
    R_obj=1738.0;
elseif strcmp(CenBody, 'Mars')==1
    mu=42828.314;
    R_obj=3396.19;
elseif strcmp(CenBody, 'Jupiter')==1
    mu=126712767.858;
    R_obj=71492;
elseif strcmp(CenBody, 'Saturn')==1
    mu=37940626.061;
    R_obj=60268;
elseif strcmp(CenBody, 'Titan')==1
    mu=8978.138;
    R_obj=2574.7;
elseif strcmp(CenBody, 'Uranus')==1
    mu=5794549.007;
    R_obj=25559;
elseif strcmp(CenBody, 'Neptune')==1
    mu=6836534.064;
    R_obj=24764;
end