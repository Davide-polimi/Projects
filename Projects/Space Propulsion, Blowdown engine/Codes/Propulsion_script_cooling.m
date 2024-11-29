%% SPACE PROPULSION PROJECT with COOLING

%%%%%%%%%%%%%%%%%%%%%%%% HEETS TRANSFER PROJECT %%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

% Data from assignment
P_c = 50e5;                             % Chamber pressure [Pa]
T_in = 1000;                            % Thrust [N]
P_c_min = 20e5;                         % Minimum Chamber pressure [Pa]
r_avail = 0.5;                          % Radius [m]
h_avail = 2;                            % Height [m]
% Oxidizer & Fuel Data:
temp_RP1 = 298.14;                      % RP-1 Storage Temperature [K]
ent_RP1 = -5430;                        % RP-1 Storage Enthalpy [cal/mol]   
rho_fuel = 580;                         % RP-1 Density [kg/m^3]
temp_LOX = 90.14;                       % LOX Storage Temperature [K]
ent_LOX = -3032;                        % LOX Storage Enthalpy [cal/mol]
rho_ox = 1140;                          % LOX density [kg/m^3]

% Constants
g0 = 9.81;                              % [m/s^2]
R = 8.314;                              % [J/mol*K]
gamma_He = 5/3;                         % Helium gamma

% Sizing Inputs:
cooling = "NO";
if cooling == "SI"
    OF = 2;                             % O/F ratio for cooling [-]
else
    OF = 2.2778;
end

% Constraints/Assumptions
eps_exit = 160;                         % Nozzle Expansion Ratio [-]
eps_c = 9;                              % Contraction ratio [-]
L_star = 1.08;                          % Characteristic CC length [m]
alpha = deg2rad(15);                    % Aperture angle reference conical nozzle [rad]
beta = deg2rad(45);                     % Closure angle reference conical nozzle [rad]
M_c = 0.066;                            % Chamber Mach number with previous analysis with CEA [-]
B_opt = 3.1;                            % Study on B

%% Engine Sizing:
[engineSize,engineSizeCEA] = engineSizing(P_c,T_in,OF,eps_c,eps_exit,L_star);

mdot = engineSize.mass.mdot;            % [kg/s]
mdot_ox = engineSize.mass.mdotLOX;      % [kg/s]
mdot_fuel = engineSize.mass.mdotRP1;    % [kg/s]

%% Cooling Sizing

T_limit_RP1 = 650;                      % Cooling limit temperature [K]

[coolingSize] = coolingSizing(engineSize,engineSizeCEA,T_limit_RP1,mdot_fuel);

%% Pressure Losses

v_pipe = 10;                            % Lines velocity [m/s]
dP_dyn_ox = 0.5*rho_ox*v_pipe^2;
dP_dyn_fuel = 0.5*rho_fuel*v_pipe^2;
dP_feed = 50e3+20e3;                    % Losses due to line (50) & valves (20) [Pa] 
dP_inj = 0.2*P_c;
if cooling == "SI"
    dP_cool = coolingSize.losses.DP_cool;
else
    dP_cool = 0;
end

%% Tanks sizing

%% BLOW DOWN OPTIMIZATION
d_Blow = 0.1;
Blow_min = 1.5;
Blow_max = 6;
ind = 0;
max_time_v = [];
Blow_vec = Blow_min:d_Blow:Blow_max;

% Blow-down
for j = Blow_vec

B = j;
ind = ind +1;

% Titanium
material.name = "Titanium-T6Al4V";
material.sigma_y.amb = 0.7612e9;        % yeld stress
material.e_w.amb = 0.9;                 % jucture weld efficiency
material.rho.amb = 4500;                % kg/m^3
material.sigma_y.crio = 1.099e9;        % yeld stress
material.e_w.crio = 0.9;                % jucture weld efficiency
material.rho.crio = 4500;               % kg/m^3



tank.P_ox = P_c+dP_inj+dP_feed+dP_dyn_ox;
tank.P_fuel = P_c+dP_inj+dP_feed+dP_dyn_fuel+dP_cool;
tank.P_max = max(tank.P_fuel,tank.P_ox);
tank.V_ratio = mdot_ox/rho_ox*rho_fuel/mdot_fuel;
tank.pipes_empty = 0.02;                % Empty space structure-tank (for pipes passage)

Gas.name= "Helium";                     % HELIUM as pressurising gas
Gas.R = 2077.3;                         % [J/ kgK]
Gas.k = gamma_He;                           

[tank] = cyl_sizing(B_opt,tank,r_avail,h_avail,engineSize,material, ...
                    Gas,rho_ox,rho_fuel,temp_LOX,temp_RP1);

%% Injectors sizing
Cd_ox = 0.65;
Cd_fuel = 0.65;
tol_inj = 5e-5;                 % Tolerance on injector diameter
tol_Cd = 0.068;                 % Tolerance on discharge coefficient

[geom_inj] = injectors(dP_inj,mdot_ox,mdot_fuel,rho_fuel,rho_ox,Cd_ox,Cd_fuel);

%% Pipes Sizing

A_pipe_fuel = mdot_fuel/(rho_fuel*v_pipe);
A_pipe_ox = mdot_ox/(rho_ox*v_pipe);
% Darcy factor losses
D_f_fuel = dP_feed*pi^2*rho_fuel/(8*mdot_fuel^2);
D_f_ox = dP_feed*pi^2*rho_ox/(8*mdot_ox^2);

%% Initialization

% Time
t_end = 3000;                                       % [s]
dt = 1;
t_simulation = 0:dt:t_end+dt;

% Coefficients

K_inj_ox = (1/(2*rho_ox))*(1/(geom_inj.A_ox_tot*Cd_ox))^2;
K_idr_ox = 8*D_f_ox/(pi^2*rho_ox);
K_dyn_ox = 1/(2*rho_ox*A_pipe_ox^2);

K_inj_fuel = (1/(2*rho_fuel))*(1/(geom_inj.A_fu_tot*Cd_fuel))^2;
K_idr_fuel = 8*D_f_fuel/(pi^2*rho_fuel);
K_dyn_fuel = 1/(2*rho_fuel *A_pipe_fuel^2);

%% Nominal Cycle

[bd_nominal] = blowdown(engineSize,engineSizeCEA,t_end,dt,K_inj_ox, ...
                        K_idr_ox,K_dyn_ox,K_inj_fuel,K_idr_fuel,K_dyn_fuel, ...
                        T_in,P_c,P_c_min,eps_exit,temp_RP1,temp_LOX,ent_RP1, ...
                        ent_LOX,rho_fuel,rho_ox,R,gamma_He,tank,cooling, ...
                        coolingSize);

bd_nominal.burning_time = bd_nominal.t_burn(find(bd_nominal.t_burn==0,1)-1);

%% Thrust Inferiore (tolleranze sull'area degli iniettori)

d_inj_fuel_inf = geom_inj.D_fu_single-tol_inj;
d_inj_ox_inf = geom_inj.D_ox_single-tol_inj;

A_inj_fuel_inf = geom_inj.N_fu*pi*d_inj_fuel_inf^2/4;
A_inj_ox_inf = geom_inj.N_ox*pi*d_inj_ox_inf^2/4;

K_inj_ox_inf = (1/(2*rho_ox))*(1/(A_inj_ox_inf*Cd_ox))^2;
K_inj_fuel_inf = (1/(2*rho_fuel))*(1/(A_inj_fuel_inf*Cd_fuel))^2;

[bd_Ainj_inf] = blowdown(engineSize,engineSizeCEA,t_end,dt,K_inj_ox_inf, ...
                        K_idr_ox,K_dyn_ox,K_inj_fuel_inf,K_idr_fuel,K_dyn_fuel, ...
                        T_in,P_c,P_c_min,eps_exit,temp_RP1,temp_LOX,ent_RP1, ...
                        ent_LOX,rho_fuel,rho_ox,R,gamma_He,tank,cooling, ...
                        coolingSize);

bd_Ainj_inf.burning_time = bd_Ainj_inf.t_burn(find(bd_Ainj_inf.t_burn==0,1)-1);

%% Thrust Superiore (tolleranze sull'area degli iniettori)

d_inj_fuel_sup = geom_inj.D_fu_single+tol_inj;
d_inj_ox_sup = geom_inj.D_ox_single+tol_inj;

A_inj_fuel_sup = geom_inj.N_fu*pi*d_inj_fuel_sup^2/4;
A_inj_ox_sup = geom_inj.N_ox*pi*d_inj_ox_sup^2/4;

K_inj_ox_sup = (1/(2*rho_ox))*(1/(A_inj_ox_sup*Cd_ox))^2;
K_inj_fuel_sup = (1/(2*rho_fuel))*(1/(A_inj_fuel_sup*Cd_fuel))^2;

[bd_Ainj_sup] = blowdown(engineSize,engineSizeCEA,t_end,dt,K_inj_ox_sup, ...
                        K_idr_ox,K_dyn_ox,K_inj_fuel_sup,K_idr_fuel,K_dyn_fuel, ...
                        T_in,P_c,P_c_min,eps_exit,temp_RP1,temp_LOX,ent_RP1, ...
                        ent_LOX,rho_fuel,rho_ox,R,gamma_He,tank,cooling, ...
                        coolingSize);

bd_Ainj_sup.burning_time = bd_Ainj_sup.t_burn(find(bd_Ainj_sup.t_burn==0,1)-1);

%% Thrust Inferiore (tolleranze sul Cd)

Cd_ox_inf = Cd_ox-tol_Cd;
Cd_fuel_inf = Cd_fuel-tol_Cd;

K_inj_ox_inf = (1/(2*rho_ox))*(1/(geom_inj.A_ox_tot*Cd_ox_inf))^2;
K_inj_fuel_inf = (1/(2*rho_fuel))*(1/(geom_inj.A_fu_tot*Cd_fuel_inf))^2;

[bd_Cd_inf] = blowdown(engineSize,engineSizeCEA,t_end,dt,K_inj_ox_inf, ...
                        K_idr_ox,K_dyn_ox,K_inj_fuel_inf,K_idr_fuel,K_dyn_fuel, ...
                        T_in,P_c,P_c_min,eps_exit,temp_RP1,temp_LOX,ent_RP1, ...
                        ent_LOX,rho_fuel,rho_ox,R,gamma_He,tank,cooling, ...
                        coolingSize);

bd_Cd_inf.burning_time = bd_Cd_inf.t_burn(find(bd_Cd_inf.t_burn==0,1)-1);

%% Thrust Superiore (tolleranze sul Cd)

Cd_ox_sup = Cd_ox+tol_Cd;
Cd_fuel_sup = Cd_fuel+tol_Cd;

K_inj_ox_sup = (1/(2*rho_ox))*(1/(geom_inj.A_ox_tot*Cd_ox_sup))^2;
K_inj_fuel_sup = (1/(2*rho_fuel))*(1/(geom_inj.A_fu_tot*Cd_fuel_sup))^2;

[bd_Cd_sup] = blowdown(engineSize,engineSizeCEA,t_end,dt,K_inj_ox_sup, ...
                        K_idr_ox,K_dyn_ox,K_inj_fuel_sup,K_idr_fuel,K_dyn_fuel, ...
                        T_in,P_c,P_c_min,eps_exit,temp_RP1,temp_LOX,ent_RP1, ...
                        ent_LOX,rho_fuel,rho_ox,R,gamma_He,tank,cooling, ...
                        coolingSize);

bd_Cd_sup.burning_time = bd_Cd_sup.t_burn(find(bd_Cd_sup.t_burn==0,1)-1);

%% Plot 

close all;
set(groot, 'defaultTextInterpreter', 'latex')
set(groot, 'defaultAxesXMinorGrid','on','defaultAxesXMinorGridMode','manual');
set(groot, 'defaultAxesYMinorGrid','on','defaultAxesYMinorGridMode','manual');
set(groot, 'defaultLegendInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter', 'latex')
set(groot, 'defaultAxesFontWeight', 'bold')
set(groot, 'defaultFigurePosition', [470, 360, 700, 430])
set(groot, 'defaultFigureColormap', turbo(256));
set(groot, 'defaultAxesFontName', 'Palatino Linotype', 'defaultTextFontName', 'Palatino Linotype');
set(groot, 'defaultSurfaceEdgeAlpha', 0.3);
set(groot, 'defaultLineLineWidth', 1.6);
set(groot, 'defaultFigureColor', [1; 1; 1]);
set(groot, 'defaultAxesColor', 'none');
set(groot, 'defaultAxesFontSize', 20);
 
figure(1)
plot(0:bd_nominal.burning_time,bd_nominal.T(1:bd_nominal.burning_time+1))
grid on
hold on
plot(0:bd_Ainj_inf.burning_time,bd_Ainj_inf.T(1:bd_Ainj_inf.burning_time+1))
plot(0:bd_Ainj_sup.burning_time,bd_Ainj_sup.T(1:bd_Ainj_sup.burning_time+1))
plot(0:bd_Cd_inf.burning_time,bd_Cd_inf.T(1:bd_Cd_inf.burning_time+1))
plot(0:bd_Cd_sup.burning_time,bd_Cd_sup.T(1:bd_Cd_sup.burning_time+1))
xlabel('$t$ $[s]$')
ylabel('$T$ $[N]$')
legend('$Nominal$ $Value$','$Maximum$ $value$ $from$ $Injectors$ $tolerances$', ...
        '$Minimum$ $value$ $from$ $Injectors$ $tolerances$', ...
        '$Maximum$ $value$ $from$ $C_d$ $tolerances$', ...
        '$Minimum$ $value$ $from$ $C_d$ $tolerances$', ...
        'Location','best')

figure(2)
plot(0:bd_nominal.burning_time,bd_nominal.P_cc(1:bd_nominal.burning_time+1))
grid on
hold on
plot(0:bd_Ainj_inf.burning_time,bd_Ainj_inf.P_cc(1:bd_Ainj_inf.burning_time+1))
plot(0:bd_Ainj_sup.burning_time,bd_Ainj_sup.P_cc(1:bd_Ainj_sup.burning_time+1))
plot(0:bd_Cd_inf.burning_time,bd_Cd_inf.P_cc(1:bd_Cd_inf.burning_time+1))
plot(0:bd_Cd_sup.burning_time,bd_Cd_sup.P_cc(1:bd_Cd_sup.burning_time+1))
plot(t_simulation,20*10^5*ones(1,length(t_simulation)),'r')
xlabel('$t$ $[s]$')
ylabel('$P$ $[Pa]$')
ylim([15*10^5 55*10^5])
legend('$Nominal$ $Value$','$Maximum$ $value$ $from$ $Injectors$ $tolerances$', ...
        '$Minimum$ $value$ $from$ $Injectors$ $tolerances$', ...
        '$Maximum$ $value$ $from$ $C_d$ $tolerances$', ...
        '$Minimum$ $value$ $from$ $C_d$ $tolerances$', ...
        '$Minimum$ $Chamber$ $Pressure$','Location','best')

figure(3)
plot(0:bd_nominal.burning_time,bd_nominal.O_F(1:bd_nominal.burning_time+1))
grid on
hold on
plot(0:bd_Ainj_inf.burning_time,bd_Ainj_inf.O_F(1:bd_Ainj_inf.burning_time+1))
plot(0:bd_Ainj_sup.burning_time,bd_Ainj_sup.O_F(1:bd_Ainj_sup.burning_time+1))
plot(0:bd_Cd_inf.burning_time,bd_Cd_inf.O_F(1:bd_Cd_inf.burning_time+1))
plot(0:bd_Cd_sup.burning_time,bd_Cd_sup.O_F(1:bd_Cd_sup.burning_time+1))
xlabel('$t$ $[s]$')
ylabel('$O/F$ $[-]$')
legend('$Nominal$ $Value$','$Maximum$ $value$ $from$ $Injectors$ $tolerances$', ...
        '$Minimum$ $value$ $from$ $Injectors$ $tolerances$', ...
        '$Maximum$ $value$ $from$ $C_d$ $tolerances$', ...
        '$Minimum$ $value$ $from$ $C_d$ $tolerances$', ...
        'Location','best')

figure(4)
plot(0:bd_nominal.burning_time,bd_nominal.mdot_v(1:bd_nominal.burning_time+1))
grid on
hold on
plot(0:bd_Ainj_inf.burning_time,bd_Ainj_inf.mdot_v(1:bd_Ainj_inf.burning_time+1))
plot(0:bd_Ainj_sup.burning_time,bd_Ainj_sup.mdot_v(1:bd_Ainj_sup.burning_time+1))
plot(0:bd_Cd_inf.burning_time,bd_Cd_inf.mdot_v(1:bd_Cd_inf.burning_time+1))
plot(0:bd_Cd_sup.burning_time,bd_Cd_sup.mdot_v(1:bd_Cd_sup.burning_time+1))
xlabel('$t$ $[s]$')
ylabel('$\dot{m}$ $[kg/s]$')
legend('$Nominal$ $Value$','$Maximum$ $value$ $from$ $Injectors$ $tolerances$', ...
        '$Minimum$ $value$ $from$ $Injectors$ $tolerances$', ...
        '$Maximum$ $value$ $from$ $C_d$ $tolerances$', ...
        '$Minimum$ $value$ $from$ $C_d$ $tolerances$', ...
        'Location','best')

if cooling == "SI"
    figure(5)
    plot(0:bd_nominal.burning_time,bd_nominal.cooling.T_end_cooling(1:bd_nominal.burning_time+1))
    grid on
    hold on
    plot(0:bd_Ainj_inf.burning_time,bd_Ainj_inf.cooling.T_end_cooling(1:bd_Ainj_inf.burning_time+1))
    plot(0:bd_Ainj_sup.burning_time,bd_Ainj_sup.cooling.T_end_cooling(1:bd_Ainj_sup.burning_time+1))
    plot(0:bd_Cd_inf.burning_time,bd_Cd_inf.cooling.T_end_cooling(1:bd_Cd_inf.burning_time+1))
    plot(0:bd_Cd_sup.burning_time,bd_Cd_sup.cooling.T_end_cooling(1:bd_Cd_sup.burning_time+1))
    plot(t_simulation,650*ones(1,length(t_simulation)),'r')
    xlabel('$t$ $[s]$')
    ylabel('$T_RP1$ $[K]$')
    legend('$Nominal$ $Value$','$Maximum$ $value$ $from$ $Injectors$ $tolerances$', ...
        '$Minimum$ $value$ $from$ $Injectors$ $tolerances$', ...
        '$Maximum$ $value$ $from$ $C_d$ $tolerances$', ...
        '$Minimum$ $value$ $from$ $C_d$ $tolerances$', ...
        'Location','best')
end


