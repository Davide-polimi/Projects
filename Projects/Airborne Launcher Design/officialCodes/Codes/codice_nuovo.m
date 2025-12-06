clc;
close all;
clear;

%% DELTAV estimation
% Data
R_E = 6378e3; % m
h_orbit = 400e3; % m
g0 = 9.81; % m/s^2
h_launch = 12e3; % m
g = g0/(1 + h_launch/R_E)^2; % m/s^2

v_req = sqrt(g * R_E^2/(R_E+h_orbit)); % circular orbit approximation

% overexstimation - consider New Zealand
lambda = deg2rad(40); % rad
v_E = 7.2921159e-5 * R_E; % m/s
v_site_loss = v_E * cos(lambda);

% gravity loss

v_grav_loss = sqrt(2 * g * R_E^2*(1/(R_E + h_launch) - 1/(R_E + h_orbit)));

% aircraft velocity
T = -55 + 273.15; % K
gamma = 1.4; % -
R = 287.02; % J/(Kg * K)
M = 0.82; % -
a = sqrt(gamma * R * T); % m/s
v_aircraft = M * a;

% Consider a percentage because deltaT_ignition = 5 s

v_budget = v_site_loss + v_grav_loss + v_req - v_aircraft;


%% INPUT DATA

OF = 2.3; % Oxidizer to fuel ratio
rho_ox = 1142; % kg/m^3
rho_fu = 810; % kg/m^3
t = 7e-2; % Tank thickness [m]
epsilon_s = 0.06; % Structural mass index
L_fair = 3.63; % Fairing length [m]
M_fairing = 145; % Fairing mass [kg]
M_pay = 250; % kg in SSO


%% Dati dei singoli stadi (modifica per adattarsi a due stadi)
stages = {
    struct('Thrust', 327e3, 'Isp', 309, 'tb', 190, 'D', 1.8), ... % First stage 
    struct('Thrust', 22.24e3, 'Isp', 328, 'tb', 400, 'D', 1.5)    % Second stage
};

%% Funzione per calcolare le proprietà di ogni stadio in base alla massa
calc_stage = @(Thrust, Isp, tb) struct( ...
    'm_flow', Thrust / (Isp * g0), ... % Propellant mass flow rate (ASSUMED CONSTANT)
    'M_prop', Thrust / (Isp * g0) * tb, ... % Propellant mass
    'M_fu', @(M_prop) OF / (1 + OF) * M_prop, ... % Fuel mass
    'M_ox', @(M_prop) 1 / (1 + OF) * M_prop ... % Oxidizer mass
);

%% Calcolo delle proprietà dei singoli stadi e delle densità medie
stage_props = cellfun(@(s) calc_stage(s.Thrust, s.Isp, s.tb), stages, 'UniformOutput', false);

for i = 1:length(stage_props)
    % Masse di ossidante e carburante per ciascun stadio
    M_prop = stage_props{i}.M_prop;
    M_fu = stage_props{i}.M_fu(M_prop);
    M_ox = stage_props{i}.M_ox(M_prop);

    % Densità media per serbatoio
    stage_props{i}.rho_av = (M_fu + M_ox) / (M_fu / rho_fu + M_ox / rho_ox);

    % Volume del serbatoio (in base a una stima razionale)
    stage_props{i}.Volume = M_prop / stage_props{i}.rho_av;
end

I_sp_av = (stage_props{1}.M_prop*stages{1}.Isp + stages{2}.Isp * stage_props{2}.m_flow * stages{1}.tb)/(stage_props{1}.M_prop + stage_props{2}.m_flow * stages{1}.tb);

% compute M0/Mf target
R_targ = exp(v_budget/(I_sp_av * g0)); 
%% Calcolo dei diametri e delle lunghezze dei serbatoi

eps_guess = (0.05:0.001:0.15);
tolerance = 1e-1; % Tolerance
max_iterations = 100; % MAximum number of iterations
rho_launcher = 610; % kg/m^2
counter = 0;
converged = false;

while ~converged && counter < max_iterations
 counter = counter + 1;
  for i = 1:length(eps_guess)
    % Massa strutturale per ciascun stadio
    stage_props{1}.M_structure = stage_props{1}.M_prop * eps_guess(i) / (1 - eps_guess(i));
    stage_props{2}.M_structure = stage_props{2}.M_prop * eps_guess(i) / (1 - eps_guess(i));
    
   
    M0_noPay = sum([stage_props{1}.M_prop, stage_props{2}.M_prop]) + ...
    sum([stage_props{1}.M_structure, stage_props{2}.M_structure]);

     % Verifica della convergenza
   R = (M0_noPay + M_pay)/M_pay;
    if abs(R - R_targ) < tolerance
        epsilon = eps_guess(i);
        converged = true;
        break;
    end
    
  end
end


