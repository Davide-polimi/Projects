clc
close all
clear

%% DELTA V estimation

% Data
R_E = 6378e3; % m
h = 400e3; % m
g0 = 9.81; % m/s^2
z = 12e3; % m
g = g0/(1 + z/R_E)^2; % m/s^2

v_req = sqrt(g * R_E^2/(R_E+h)); % circular orbit approximation

% overexstimation - consider New Zealand
lambda = deg2rad(40); % rad
v_E = 7.2921159e-5 * R_E; % m/s
v_site_loss = v_E * cos(lambda);

% gravity loss

v_grav_loss = sqrt(2 * g * R_E^2*(1/(R_E + z) - 1/(R_E + h)));

% aircraft velocity
T = -55 + 273.15; % K
gamma = 1.4; % -
R = 287.02; % J/(Kg * K)
M = 0.82; % -
a = sqrt(gamma * R * T); % m/s
v_aircraft = M * a;

% Consider a percentage because deltaT_ignition = 5 s

v_budget = v_site_loss + v_grav_loss + v_req - v_aircraft;
I_sp = 300;

MR = exp(v_budget/(I_sp * g0));

%% Total mass estimation
% LauncherONE DATA % to be updated............ 
Th1 = 238.6e3; % N
Th2 = 22.24e3; % N
Is1 = 309; % s
Is2 = 328; % s
tb1 = 190; % s
tb2 = 400; % s
D1 = 1.7; % m
D2 = 1.7; % m
D_fair = 1.7; % m
L_fair = 4.57; % m

% Efficiency of Engine, compute the ideal specific impulses from CEA code

Is_id1 = 339; % s
Is_id2 = 364; % s

eff_1 = Is1/Is_id1;
eff_2 = Is2/Is_id2;

% Mass flow rates
m_flow_1 = Th1/(Is1*g0);
m_flow_2 = Th2/(Is2*g0);

% Propellant masses

M_p1 = m_flow_1 * tb1;
M_p2 = m_flow_2 * tb2;

rho_ox = 1142; % kg/m^3
rho_fu = 810; % kg/m^3
OF = 2.3;

M_fu_1 = OF/(1 + OF) * M_p1; 
M_fu_2 = OF/(1 + OF) * M_p2;

M_ox_1 = 1/(1 + OF) * M_p1;
M_ox_2 = 1/(1 + OF) * M_p2;

% Average density
rho_av_1 =(M_fu_1+M_ox_1)/(M_fu_1/rho_fu + M_ox_1/rho_ox); 
rho_av_2 =(M_fu_2+M_ox_2)/(M_fu_2/rho_fu + M_ox_2/rho_ox); 

% Volumes
Vol_1 = M_p1/rho_av_1;
Vol_2 = M_p2/rho_av_2;

% suppose the thickness of the tanks

t = 7e-2; % m

L_tank1 = Vol_1*4/(pi*((D1-2*t)^2));
L_tank2 = Vol_2*4/(pi*((D2-2*t)^2));

L_tot = 21.336; % m

L_other = L_tot - (L_tank2 + L_tank1 + L_fair);

L_per = L_other/L_tot; % percentage of length occupied by CC + void

% Stage lengths

L_s1 = L_tank1/(1 - L_per);
L_s2 = L_tank2/(1 - L_per);

epsilon_s = 0.06;

M_s1 = M_p1 * epsilon_s/(1-epsilon_s);
M_s2 = M_p2 * epsilon_s/(1-epsilon_s);

% M_S = 10; % kg/m^2
M_fairing = 145; % kg
M0 = M_p1 + M_s1 + M_p2 + M_s2 + M_fairing;

%% Eventuale terzo stadio
% Parametri del terzo stadio
Th3 = 0; % N
Is3 = 0; % s
tb3 = 0; % s
D3 = 0; % m

% Efficienza del terzo stadio
Is_id3 = 0; % s
eff_3 = Is3 / Is_id3;

% Flusso di massa e massa di propellente
m_flow_3 = Th3 / (Is3 * g0);
M_p3 = m_flow_3 * tb3;

% Masse di combustibile e ossidante
M_fu_3 = OF / (1 + OF) * M_p3;
M_ox_3 = 1 / (1 + OF) * M_p3;

% Densità media e volume
rho_av_3 = (M_fu_3 + M_ox_3) / (M_fu_3 / rho_fu + M_ox_3 / rho_ox);
Vol_3 = M_p3 / rho_av_3;

% Lunghezza del serbatoio
L_tank3 = Vol_3 * 4 / (pi * ((D3 - 2 * t)^2));

% Lunghezza e percentuale altre componenti
L_tot = 21.336; % m
L_other = L_tot - (L_tank1 + L_tank2 + L_tank3 + L_fair);
L_per = L_other / L_tot;

% Massa strutturale del terzo stadio
M_s3 = M_p3 * epsilon_s / (1 - epsilon_s);

% Massa totale iniziale
M0 = M_p1 + M_s1 + M_p2 + M_s2 + M_p3 + M_s3 + M_fairing;



