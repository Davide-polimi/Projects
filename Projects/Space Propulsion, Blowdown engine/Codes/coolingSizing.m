function [coolingSize] = coolingSizing(engineSize, engineSizeCEA, T_limit_RP1, mdot_RP1)

% OUTPUTS: temperature di RP1 in:
    % - fine della camera di combustione
    % - gola
    % - 

D_c = engineSize.diameters.Dchamber;
D_t = engineSize.diameters.Dthroat;
D_e = engineSize.diameters.Dexit;
L_c = engineSize.lengths.Lchamber;
L_conv = engineSize.lengths.Lconv;
L_div = engineSize.lengths.Ldiv;
Tc = engineSizeCEA.output.froz.temperature(1);

% Oxidizer & Fuel Data:
temp_RP1 = 298;                         % RP-1 Storage Temperature [K]
ent_RP1 = -5430;                        % RP-1 Storage Enthalpy [cal/mol]   
temp_LOX = 90;                          % LOX Storage Temperature [K]
ent_LOX = -3032;                        % LOX Storage Enthalpy [cal/mol]

%%%%%%%%%%%%%%%%%%%%%%%% COOLING PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RP-1 Data:
rho_RP1 = 799.67;                       % RP-1 Density [kg/m^3]
mu_RP1 = 1.550*10^-3;                   % RP-1 Cinematic Viscosity [Pa*s]   from: 1-s2.0-S0016236118314418-main.pdf
cond_RP1 = 0.11324;                     % RP-1 Thermal Conductivity [W/mK]  from: 1-s2.0-S0016236118314418-main.pdf
cp_RP1 = 1.88*1000;                     % RP-1 Specific Heat [J/kg*K]
Pr_RP1 = cp_RP1*mu_RP1/cond_RP1;        % RP-1 Prandtl Number

% Wall Data (Incolonel X):              Ref: pg. 109 libro Huzel & Huang:
Tw = 973;                               % Wall Limit Temperature [K]         da: https://asm.matweb.com/search/SpecificMaterial.asp?bassnum=NINC35
cond_wall = 12;                         % Wall Thermal Conductivity [W/mK]   da: https://asm.matweb.com/search/SpecificMaterial.asp?bassnum=NINC35

% Thermal Barrier Coating (Y2O3):       Ref: https://www.makeitfrom.com/material-properties/Yttria-Yttrium-Oxide-Y2O3
T_TBC = 2000;                           % TBC Limit Temperature [K]         
cond_TBC = 2;                           % TBC Thermal Conductivity [W/mK]  

% Bond properties:                      Ref: https://www.researchgate.net/publication/281579071_Stress_evolution_in_thermal_Barrier_coatings_for_rocket_engine_applications
T_bond = 1300;                          % Bond Limit Temperature [K]  
cond_bond = 17;                         % Bond Thermal Conductivity [W/mK]  

%%%%%%%%%%%%%%%%%%%%%%%%%%%% CEA OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A: Combustion Chamber 
% B: Nozzle Convergent
% C: Nozzle Divergent

% Mach Number [-]:
M_A = engineSizeCEA.output.froz.mach(1);
M_B = engineSizeCEA.output.froz.mach(2);
M_C = engineSizeCEA.output.froz.mach(3);
% Sonic Velocity [m/s]:
a_A = engineSizeCEA.output.froz.sonvel(1);
a_B = engineSizeCEA.output.froz.sonvel(2);
a_C = engineSizeCEA.output.froz.sonvel(3);
% Velocity [m/s]:
v_A = M_A*a_A;
v_B = M_B*a_B;
v_C = M_C*a_C;
% Specific Heat [KJ/kg*K]:
cp_A = engineSizeCEA.output.froz.cp_tran.froz(1);
cp_B = engineSizeCEA.output.froz.cp_tran.froz(2);
cp_C = engineSizeCEA.output.froz.cp_tran.froz(3);
% Gamma [-]:
k_A = engineSizeCEA.output.froz.gamma(1);
k_B = engineSizeCEA.output.froz.gamma(2);
k_C = engineSizeCEA.output.froz.gamma(3);
% Cinematic Viscosity [Pa*s]:
mu_A = engineSizeCEA.output.froz.viscosity(1)*10^-6;
mu_B = engineSizeCEA.output.froz.viscosity(2)*10^-6;
mu_C = engineSizeCEA.output.froz.viscosity(3)*10^-6;
% Density [kg/m^3]:
rho_A = engineSizeCEA.output.froz.density(1);
rho_B = engineSizeCEA.output.froz.density(2);
rho_C = engineSizeCEA.output.froz.density(3);
% Prandtl Number:
Pr_A = engineSizeCEA.output.froz.prandtl.froz(1);
Pr_B = engineSizeCEA.output.froz.prandtl.froz(2);
Pr_C = engineSizeCEA.output.froz.prandtl.froz(3);
% Conductivity [W/mK]:
cond_A = engineSizeCEA.output.froz.conduct.froz(1);
cond_B = engineSizeCEA.output.froz.conduct.froz(2);
cond_C = engineSizeCEA.output.froz.conduct.froz(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Heat Transfer Process %%%%%%%%%%%%%%%%%%%%%%%%

% Convective Heat Transfer Coefficient (gas side):
Re_A = rho_A*v_A*D_c/mu_A;
Nu_A = 0.0265* Re_A^0.8 * Pr_A^0.3;
Re_B = rho_B*v_B*D_t/mu_B;
Nu_B = 0.0265* Re_B^0.8 * Pr_B^0.3;
Re_C = rho_C*v_C*D_e/mu_C;
Nu_C = 0.0265* Re_C^0.8 * Pr_C^0.3;

h_gA = Nu_A*cond_A/D_c;                 % Convective HT Coefficient (A) [W/m^2K]
h_gB = Nu_B*cond_B/D_t;                 % Convective HT Coefficient (B) [W/m^2K]
h_gC = Nu_C*cond_C/D_e;                 % Convective HT Coefficient (C) [W/m^2K]

% Desired Wall Temperature (coolant side) [K]:
Twc_des = T_limit_RP1;              

% Adiabatic Wall Temperature [K]:
R_A = (1+ Pr_A^(1/3) * 0.5*(k_A - 1)*M_A^2)/(1 + 0.5*(k_A - 1)*M_A^2);
R_B = (1+ Pr_B^(1/3) * 0.5*(k_B - 1)*M_B^2)/(1 + 0.5*(k_B - 1)*M_B^2);
R_C = (1+ Pr_C^(1/3) * 0.5*(k_C - 1)*M_C^2)/(1 + 0.5*(k_C - 1)*M_C^2);

Taw_A = R_A * Tc;                       % Adiabatic Wall Temperature (A) [K]
Taw_B = R_B * Tc;                       % Adiabatic Wall Temperature (B) [K]
Taw_C = R_C * Tc;                       % Adiabatic Wall Temperature (C) [K]

% Heat Fluxes per surface unity [W/m^2]:
qA = h_gA*(Taw_A - T_TBC);              % Heat Flux per surface unity (A) [W/m^2]
qB = h_gB*(Taw_B - T_TBC);              % Heat Flux per surface unity (B) [W/m^2]
qC = h_gC*(Taw_C - T_TBC);              % Heat Flux per surface unity (C) [W/m^2]

% Side Areas [m^2]:
Aside_cc = pi*D_c*L_c;                                                      % Side Area (A) [m^2]
Aside_conv = 0.5*pi*(D_t + D_c)*sqrt(L_conv^2 + (0.5*D_c - 0.5*D_t)^2);     % Side Area (B) [m^2]
  

% Thickness of the different Layers (fixing the Temperatures)

tTBC = cond_TBC/qA * (T_TBC - T_bond);      % TBC Thickness [m]
t_bond = cond_bond/qA * (T_bond - Tw);      % Bond Thickness [m]
twall_A = cond_wall/qA * (Tw - Twc_des);    % Wall Thickness (A) [m]
twall_B = twall_A;                          % Wall Thickness (B) [m]


%%%%%% RADIATIVE COOLING IN THE NOZZLE %%%%%%%%%
% Material: 
% Reference:
% Data:
k_boltzmann = 5.67*10^-8;
emissivity = 0.42; 
Tw_rad = T_bond - qC*t_bond/cond_bond;
q_rad_c = emissivity*k_boltzmann*(Tw_rad^4);


%%%%%%% ITERATIVE PROCESS TO SIZE THE TUBES FOR THE COOLING JACKET %%%%%%%%
T_co_A = temp_RP1;
h_co_des = qA/(Twc_des - T_co_A);

% calcoliamo ora h_co a partire dal numero di Nusselt facendo variare iterativamente 
% il diametro dei fori fino a quando la condizione h_co_Nu > h_co_des è
% soddisfatta:
perimetro_chamber = pi*D_c; 
Nguess = 10;
Dhydr = perimetro_chamber/Nguess; 
h_co_Nu = cond_RP1/Dhydr * 0.243 * Pr_RP1^0.4 * (4*mdot_RP1 / (pi^2 * D_c * mu_RP1))^0.8;
while h_co_Nu < h_co_des
     h_co_Nu = cond_RP1/Dhydr * 0.243 * Pr_RP1^0.4 * (4*mdot_RP1 / (pi^2 * D_c * mu_RP1))^0.8;
     Dhydr = Dhydr - 1e-8;
end
Dtube = Dhydr;                  % Tube Diameter [m]
tTube = 0.0005;                 % Tube Thickness [m]
Nactual = ceil((perimetro_chamber + 2*tTube)/(Dhydr + 2*tTube));    % Number of Tubes

%%%%%%%%%%%%%%%%%%% COOLING PROCESS --> Heating of RP1 %%%%%%%%%%%%%%%%%%%%%
Tinl_RP1 = temp_RP1;

% A (combustion chamber):
H_A = 1/(1/h_gA + twall_A/cond_wall + 1/h_co_Nu + tTBC/cond_TBC + t_bond/cond_bond);  % equivalent resistance A
DeltaT_A = Aside_cc*H_A*(Taw_A - Tinl_RP1)/(cp_RP1*mdot_RP1);
Tfin_A = Tinl_RP1 + DeltaT_A;           % Final Temperature of RP-1 after the Combustion Chamber [K]

% B (converging Part):
% find h_co desired for B-trait:
Tco_B = Tfin_A;
h_co_desB = qB/(Twc_des - Tco_B);

H_B = 1/(1/h_gB + twall_B/cond_wall + 1/h_co_desB + tTBC/cond_TBC + t_bond/cond_bond);  % equivalent resistance B
DeltaT_B = (Aside_conv*H_B*(Taw_B - Tco_B))/(cp_RP1*mdot_RP1);
Tfin_B = Tco_B + DeltaT_B;              % Final Temperature of RP-1 after the Converging Part [K]


%%%%%%%%%%%%%%%%%%%%%%%% CHECK COOLING COMPLIANCY %%%%%%%%%%%%%%%%%%%%%%%%%
% Final Temperature of RP-1 at the end of the Cooling Process [K]:
Tend_cooling = Tfin_B;
% Limit Temperature of RP-1 [K]:
Tmargin = 650; 
% Check Regenerative Cooling:
if Tend_cooling > Tmargin
    disp('Fuel overheated')
    if qC < q_rad_c
        disp('Nozzle cooled')
    else
        disp('Nozzle is melting')
    end
else
    disp('Cooling compliant')
    if qC < q_rad_c
        disp('Nozzle cooled')
    else
        disp('Nozzle is melting')
    end
end

coolingSize.constants.beta = 0.4;
coolingSize.constants.A_tube = pi*Dtube^2/4;
coolingSize.constants.C_D_cooling = 0.61;
Delta_P_cool = ((1-coolingSize.constants.beta^4)*mdot_RP1^2)/(2*rho_RP1*(coolingSize.constants.C_D_cooling*coolingSize.constants.A_tube)^2);

% OUTPUTS:
coolingSize.temperatures.T_endChamber = Tfin_A;
coolingSize.temperatures.T_throat = Tfin_B;
coolingSize.temperatures.T_RP1_end = Tend_cooling;
coolingSize.thickness.wall = twall_A;
coolingSize.thickness.bond = t_bond;
coolingSize.thickness.TBC = tTBC;
coolingSize.tubes.diameter = Dtube;
coolingSize.tubes.thickness = tTube;
coolingSize.tubes.Ntubes = Nactual;
coolingSize.losses.DP_cool = Delta_P_cool;
coolingSize.radiative.qC = qC;
coolingSize.radiative.qC_rad = q_rad_c;
coolingSize.radiative.Tw_rad = Tw_rad;