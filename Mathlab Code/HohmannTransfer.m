function dV=HohmannTransfer(hp1,ha1,hp2,ha2)
% Computes the dV required to transfer between two orbits given the orbits'
% radius of perigee & apogee.
mu=398600.48;
re=6378.14;
% alt to radius
rp1=hp1+re;
ra1=ha1+re;
rp2=hp2+re;
ra2=ha2+re;
%semimajor axis
a1=(rp1+ra1)/2;
a2=(rp2+ra2)/2;
at=(ra1+rp2)/2;
%orbit energy
nrg1=-mu/(2*a1);
nrgt=-mu/(2*at);
nrg2=-mu/(2*a2);
%velocities at burn points
va1=sqrt(2*(nrg1+mu/ra1));
vpt=sqrt(2*(nrgt+mu/ra1));
vat=sqrt(2*(nrgt+mu/rp2));
vp2=sqrt(2*(nrg2+mu/rp2));
%dVs
dV1=abs(vpt-va1);
dV2=abs(vp2-vat);
dV=dV1+dV2;