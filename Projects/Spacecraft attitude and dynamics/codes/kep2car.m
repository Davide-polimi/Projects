function [r,v] = kep2car(a, e, i, OM, om, th, mu)
% kep2car.m - Conversion from Keplerian elements to Cartesian coordinates
%
% PROTOTYPE:
% [r, v] = kep2car(a, e, i, OM, om, th, mu)
%
% DESCRIPTION:
% Conversion from Keplerian elements to Cartesian coordinates. Angles in
% radians.
%
% INPUT:
% a [1x1] Semi-major axis [km]
% e [1x1] Eccentricity [-]
% i [1x1] Inclination [rad]
% OM [1x1] RAAN [rad]
% om [1x1] Pericentre anomaly [rad]
% th [1x1] True anomaly [rad]
% mu [1x1] Gravitational parameter [km^3/s^2]
%
% OUTPUT:
% r [3x1] Position vector [km]
% v [3x1] Velocity vector [km/s]

%% 0
R3OM=[cos(OM) sin(OM) 0;
    -sin(OM) cos(OM) 0;
    0 0 1];

R1i=[1 0 0;
    0 cos(i) sin(i);
    0 -sin(i) cos(i)];

R3om=[cos(om) sin(om) 0;
    -sin(om) cos(om) 0;
    0 0 1];

%% 1
p=a*(1-e^2);
rmod=p/(1+e*cos(th));

%% 2
rpf=rmod*[cos(th);sin(th);0];
vpf=sqrt(mu/p)*[-sin(th);e+cos(th);0];

%% 3
T=(R3om*R1i*R3OM)';
r=T*rpf;
v=T*vpf;


