format short 
format compact

r_e=6378.14; %radius of Earth
mu_e=398600.44; %mu of Earth
r_geo=35786+r_e; %radius at geo

t_t=23*60*60; %transfer orbit period of 18 hours
isp_tug=350; %Isp of GEO tug's main thrusters, RP-1, UDMH, or better
mass_msv=4e+3; %Mass of Manned Servicing Vehicle, NASA's SEV is 3000 kg + 1000 kg payload
target_m_max=1e+4; %Maximum theoritical mass of target debris

%transfer orbit properties
t_t=18*60*60; %transfer orbit period of 18 hours
a_t=aperiod(t_t,mu_e); %semimajor axis
rp_t=a_t-r_geo; %radius of perigee, assumes apogee is roughly GEO

%dV to raise/lower perigee
dV_req=HohmannTransfer(r_geo,r_geo,rp_t,r_geo);

%Rendevous Clarke Platform - full load 0 reserve
mass_empty=target_m_max+mass_msv %non-fuel mass
fuel_res=PropMassReq(dV_req,mass_empty,isp_tug) %reserve fuel

%Rendevous Clarke Platform - full load w/ reserve
mass_bingo=mass_empty+fuel_res %mass without expending reserve
fuel_home=PropMassReq(dV_req,mass_bingo,isp_tug)

%Depart target - full load
mass_depart=mass_bingo+fuel_home
fuel_depart=PropMassReq(dV_req,mass_depart,isp_tug)

%Arrive target - no load
mass_arrive=mass_depart+fuel_depart-target_m_max
fuel_arrive=PropMassReq(dV_req,mass_arrive,isp_tug)

%Leave Clarke Platform - no load
mass_leave=mass_arrive+fuel_arrive;
fuel_leave=PropMassReq(dV_req,mass_leave,isp_tug);

%mass of MSV fully fuelled
mass_MSV=mass_leave+fuel_leave