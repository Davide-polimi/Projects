clear all
clc

rE = 6371008.7714; % Radius of earth
mu = 3.986004418e14;

% Design constants
structure.mDry = [0, 1819, 251]; % [kg] Dry mass
structure.m0 = [0, 25986, 2289]; % [kg] Initial mass
structure.mPayload = 250; % [kg] Payload mass

propulsion.Isp = [1, 300, 330]; % [s] Isp
propulsion.thrust = [1, 450000, 22000]; % [N] Thrust

launcher.launchVelocity = 241.9571; % [m/s]
launcher.launchPitch = 28; % [degrees]
launcher.launchHeight = 11000;
launcher.launchLatitude = 60; % [degrees]

stability.pitchRate = 3; % [degrees/second]
stability.pitchAccel = 1; % [degrees/second^2]

rocketDesign.a = [0 1.7 1.7];
rocketDesign.b = [0 1.7 1.7];
rocketDesign.x = [0 4.5 21];
rocketDesign.ln = rocketDesign.x(2);
rocketDesign.d_nose = rocketDesign.a(2);
rocketDesign.nose_type = 'C';

% Target Orbit
target.r = rE + 400000;
target.V = sqrt(mu/target.r);

% Simulation constants
Dt = 2;

[qMaxLoad, deltaV] = trajectory2D(Dt, structure, propulsion, launcher, stability, rocketDesign, target);