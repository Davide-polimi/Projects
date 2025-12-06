% Spacecraft Guidance and Navigation (2024/2025)
% Assignment # 2, Exercise 2
% Author: Davide Bellini



clearvars; close all; clc;
cspice_kclear()

plot_set(); % Plot settings

%% DATA
% Load kernels
cspice_furnsh('assignment02.tm');

addpath('sgp4');

% Constants
arcsec2rad = pi / (180*3600);
typerun    = 'u';  % user-provided inputs to SGP4 Matlab function
opsmode    = 'a';  % afspc approach ('air force space command')
whichconst =  72;  % WGS72 constants (radius, gravitational parameter)

% DATA 

% SMOS twop-line element
longstr1 = '1 36036U 09059A   24323.76060260  .00000600  00000-0  20543-3 0  9995';
longstr2 = '2 36036  98.4396 148.4689 0001262  95.1025 265.0307 14.39727995790658';

% Ground station: KOUROU
KOUROU.name = 'KOUROU';
KOUROU.coordinates.lat = 5.25144; % Latitude in degrees
KOUROU.coordinates.lon = -52.80466; % Longitude in degrees
KOUROU.coordinates.alt = -14.67; % Altitude in meters
KOUROU.type = 'Radar (monostatic)';
KOUROU.measurements.type = {'Az', 'El', 'Range (one-way)'};
KOUROU.noise.sigma_az_el = 125e-3; % Noise in Az, El [rad]
KOUROU.noise.sigma_range = 0.01; % Noise in Range [km]
KOUROU.min_elevation = 6; % Minimum elevation [deg]
KOUROU.measurement_frequency = 60; % Measurement frequency [s]
KOUROU.cost_per_pass = 30000; % Cost per pass [€]

% Ground station: TROLL
TROLL.name = 'TROLL';
TROLL.coordinates.lat = -72.011977; % Latitude in degrees
TROLL.coordinates.lon = 2.536103; % Longitude in degrees
TROLL.coordinates.alt = 1298; % Altitude in meters
TROLL.type = 'Radar (monostatic)';
TROLL.measurements.type = {'Az', 'El', 'Range (one-way)'};
TROLL.noise.sigma_az_el = 125e-3; % Noise in Az, El [rad]
TROLL.noise.sigma_range = 0.01; % Noise in Range [km]
TROLL.min_elevation = 0; % Minimum elevation [deg]
TROLL.measurement_frequency = 30; % Measurement frequency [s]
TROLL.cost_per_pass = 35000; % Cost per pass [€]

% Ground station: SVALBARD
SVALBARD.name = 'SVALBARD';
SVALBARD.coordinates.lat = 78.229772; % Latitude in degrees
SVALBARD.coordinates.lon = 15.407786; % Longitude in degrees
SVALBARD.coordinates.alt = 458; % Altitude in meters
SVALBARD.type = 'Radar (monostatic)';
SVALBARD.measurements.type = {'Az', 'El', 'Range (one-way)'};
SVALBARD.noise.sigma_az_el = 125e-3; % Noise in Az, El [rad]
SVALBARD.noise.sigma_range = 0.01; % Noise in Range [km]
SVALBARD.min_elevation = 8; % Minimum elevation [deg]
SVALBARD.measurement_frequency = 60; % Measurement frequency [s]
SVALBARD.cost_per_pass = 35000; % Cost per pass [€]

%% --------------- 2.1 Compute visibility windows -----------------------

% Initialize the satrec structure
satrec = twoline2rv(longstr1, longstr2, typerun,'e', opsmode, whichconst);

[year,mon,day,hr,min,sec] = invjday(satrec.jdsatepoch, satrec.jdsatepochf);
sat_epoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,hr,min,sec]);
sat_epoch_et = cspice_str2et(sat_epoch_str);

% Evaluate the TLE with sgp4
[satrec,rteme,vteme] = sgp4(satrec, 0.0);
% Osculating orbital elements
elts = cspice_oscelt( [rteme;vteme], sat_epoch_et, satrec.mu );

% Nutation and Precession corrections
ddpsi = -0.111690*arcsec2rad; %  [rad]
ddeps = -0.006855*arcsec2rad; %  [rad]

% Centuries from TDT 2000 January 1 00:00:00.000 
et_ref = cspice_str2et(sat_epoch_str); % Epoch we want to perform our conversion TEME to ECI
ttt = cspice_unitim(et_ref, 'ET', 'TDT')/cspice_jyear()/100;

% Transform TEME to ECI vectors
ateme = [0;0;0];
[reci, veci, aeci] = teme2eci(rteme, vteme, ateme, ttt, ddpsi, ddeps);

% Time
str_t0 = '2024-11-18T20:30:00.000';
str_tf = '2024-11-18T22:15:00.000';
et_0 = cspice_str2et(str_t0); 
et_f = cspice_str2et(str_tf);

% Initialization
y_ref = [reci;veci];
tspan_0 = [et_ref et_0];
npoints = round((et_f-et_0)/30.0) + 1;
tspan_1 = linspace(et_0, et_f, npoints);
% Options
options = odeset('RelTol',1e-13,'AbsTol',1e-20);

% Propagation from reference to starting time
[tt_0,yy_0] = ode113(@(t,y) ode_2bp(t,y,satrec.mu), tspan_0, y_ref,options);
% Propagation from starting time to final time
[tt,yy_ECI] = ode113(@(t,y) ode_2bp(t,y,satrec.mu), tspan_1, yy_0(end,1:6),options);

% Measurments
stationName = KOUROU.name;
[azimuth_Kou, elevation_Kou, ~, ~] = antennaPointing(stationName, tt(1:2:end), yy_ECI(1:2:end,:));

stationName = TROLL.name;
[azimuth_Tro, elevation_Tro, ~, ~] = antennaPointing(stationName, tt, yy_ECI);

stationName = SVALBARD.name;
[azimuth_Sva, elevation_Sva, ~, ~] = antennaPointing(stationName, tt(1:2:end), yy_ECI(1:2:end,:));


% Visibility indices
i_vis_Kou = elevation_Kou > deg2rad(KOUROU.min_elevation);
i_vis_Tro = elevation_Tro > deg2rad(TROLL.min_elevation);
i_vis_Sva = elevation_Sva > deg2rad(SVALBARD.min_elevation);

% Plot azimuth elevation-----------------------------------------------
figure
plot(azimuth_Kou(i_vis_Kou)*cspice_dpr(), elevation_Kou(i_vis_Kou)*cspice_dpr(),'*','DisplayName', 'SMOS','MarkerSize',10,'LineWidth',1.2)
hold on
plot(azimuth_Tro(i_vis_Tro)*cspice_dpr(), elevation_Tro(i_vis_Tro)*cspice_dpr(),'*','DisplayName', 'SMOS','MarkerSize',10,'LineWidth',1.2)
plot(azimuth_Sva(i_vis_Sva)*cspice_dpr(), elevation_Sva(i_vis_Sva)*cspice_dpr(),'*','DisplayName', 'SMOS','MarkerSize',10,'LineWidth',1.2)
grid on;
axis([-180,180,0, 90])
xlabel('Azimuth [deg]')
ylabel('Elevation [deg]')
legend('Kourou (60 s)','Troll (30 s)','Svalbard (60 s)')
title('Visibility Windows')


% Azimuth and elevation KOUROU-------------------------------------------
figure
yyaxis left;
plot(tt(1:2:end)/cspice_spd, elevation_Kou*cspice_dpr(), '-', 'LineWidth', 1.2);
ylabel('Elevation [deg]'); % Label per l'asse y di sinistra
hold on;
% Right y-axis
yyaxis right;
plot(tt(1:2:end)/cspice_spd, azimuth_Kou*cspice_dpr(), '-', 'LineWidth', 1.2);
ylabel('Azimuth [deg]'); % Label per l'asse y di destra
grid on;
xlabel('18-NOV-2024')
legend('Elevation','Azimuth')
title('KOUROU')
num_ticks = 4;
ax = gca;
ax.FontSize = 15;
tick_indices = round(linspace(1, length(tt), num_ticks)); 
tick_values = tt(tick_indices) / cspice_spd();
tick_labels = cell(num_ticks, 1); 
for i = 1:num_ticks
    utc_full = cspice_et2utc(tt(tick_indices(i)), 'C', 0); 
    tick_labels{i} = utc_full(12:end); 
end
xticks(tick_values);
xticklabels(tick_labels);
xtickangle(0);

% Azimuth and elevation Svalbard------------------------------------------
figure
yyaxis left;
plot(tt(1:2:end)/cspice_spd, elevation_Sva*cspice_dpr(), '-', 'LineWidth', 1.2);
ylabel('Elevation [deg]'); % Label per l'asse y di sinistra
hold on;
% Right y-axis
yyaxis right;
plot(tt(1:2:end)/cspice_spd, azimuth_Sva*cspice_dpr(), '-', 'LineWidth', 1.2);
ylabel('Azimuth [deg]'); % Label per l'asse y di destra
grid on;
xlabel('18-NOV-2024')
legend('Elevation','Azimuth')
title('SVALBARD')
num_ticks = 4;
ax = gca;
ax.FontSize = 15;
tick_indices = round(linspace(1, length(tt), num_ticks)); 
tick_values = tt(tick_indices) / cspice_spd();
tick_labels = cell(num_ticks, 1); 
for i = 1:num_ticks
    utc_full = cspice_et2utc(tt(tick_indices(i)), 'C', 0); 
    tick_labels{i} = utc_full(12:end); 
end
xticks(tick_values);
xticklabels(tick_labels);
xtickangle(0);

% Azimuth and elevation Troll-----------------------------------------
figure
yyaxis left;
plot(tt(1:end)/cspice_spd, elevation_Tro*cspice_dpr(), '-', 'LineWidth', 1.2);
ylabel('Elevation [deg]'); % Label per l'asse y di sinistra
hold on;
% Right y-axis
yyaxis right;
plot(tt(1:end)/cspice_spd, azimuth_Tro*cspice_dpr(), '-', 'LineWidth', 1.2);
ylabel('Azimuth [deg]'); % Label per l'asse y di destra
grid on;
xlabel('18-NOV-2024')
legend('Elevation','Azimuth')
title('TROLL')
num_ticks = 4;
ax = gca;
ax.FontSize = 15;
tick_indices = round(linspace(1, length(tt), num_ticks)); 
tick_values = tt(tick_indices) / cspice_spd();
tick_labels = cell(num_ticks, 1); 
for i = 1:num_ticks
    utc_full = cspice_et2utc(tt(tick_indices(i)), 'C', 0); 
    tick_labels{i} = utc_full(12:end); 
end
xticks(tick_values);
xticklabels(tick_labels);
xtickangle(0);


% VIsibility
% Visibility Time
t_vis_Tro = tt(i_vis_Tro == 1);
% Time vector for 60s sampling
tt_60 = tt(1:2:end);
t_vis_Kou = tt_60(i_vis_Kou == 1);
t_vis_Sva = tt_60(i_vis_Sva == 1);

% Initial
Vis.Kou_vis_i = cspice_et2utc(t_vis_Kou(1),'C',3);
Vis.Tro_vis_i = cspice_et2utc(t_vis_Tro(1),'C',3);
Vis.Sva_vis_i = cspice_et2utc(t_vis_Sva(1),'C',3);
% Final
Vis.Kou_vis_f = cspice_et2utc(t_vis_Kou(end),'C',3);
Vis.Tro_vis_f = cspice_et2utc(t_vis_Tro(end),'C',3);
Vis.Sva_vis_f = cspice_et2utc(t_vis_Sva(end),'C',3);
% Visibility



%% 2.2 Simulate measurements

%------------------------- KOUROU ----------------------------
% Initialization
reci_sc_Kou = zeros(3,length(t_vis_Kou));
veci_sc_Kou = zeros(3,length(t_vis_Kou));

for i = 1:length(t_vis_Kou)
    % SGP4 propagation
    tsince = (t_vis_Kou(i) - sat_epoch_et)/60.0; % minutes from TLE epoch
    [~,rteme_sc,vteme_sc] = sgp4(satrec,  tsince);
    
    % Precession
    ttt = cspice_unitim(t_vis_Kou(i), 'ET', 'TDT')/cspice_jyear()/100;
    
    % TEME to ECI conversion
    [reci_sc_Kou(:,i), veci_sc_Kou(:,i)] = ...
        teme2eci(rteme_sc, vteme_sc, [0.0;0.0;0.0],  ttt, ddpsi, ddeps);
end

stationName = KOUROU.name;
noise = KOUROU.noise;
% SGP4 without noise
[az_Kou, elev_Kou, ~, ~] = antennaPointing(stationName, t_vis_Kou, [reci_sc_Kou; veci_sc_Kou]');
% SGP4 with noise
[az_Kou_meas, elev_Kou_meas, range_Kou_meas, rangeRate_Kou_meas] = antennaPointing(stationName, t_vis_Kou, [reci_sc_Kou; veci_sc_Kou]',noise);


%----------------- TROLL ------------------------------------------------
% Initialization
reci_sc_Tro = zeros(3,length(t_vis_Tro));
veci_sc_Tro = zeros(3,length(t_vis_Tro));

for i = 1:length(t_vis_Tro)
    % SGP4 propagation
    tsince = (t_vis_Tro(i) - sat_epoch_et)/60.0; % minutes from TLE epoch
    [~,rteme_sc,vteme_sc] = sgp4(satrec,  tsince);

    % Precession
    ttt = cspice_unitim(t_vis_Tro(i), 'ET', 'TDT')/cspice_jyear()/100;

    % TEME to ECI conversion
    [reci_sc_Tro(:,i), veci_sc_Tro(:,i)] = ...
        teme2eci(rteme_sc, vteme_sc, [0.0;0.0;0.0],  ttt, ddpsi, ddeps);
end

stationName = TROLL.name;
noise = TROLL.noise;
% SGP4 without noise
[az_Tro, elev_Tro, ~, ~] = antennaPointing(stationName, t_vis_Tro,[reci_sc_Tro; veci_sc_Tro]');
% SGP4 with noise
[az_Tro_meas, elev_Tro_meas, range_Tro_meas, rangeRate_Tro_meas] = antennaPointing(stationName, t_vis_Tro,[reci_sc_Tro; veci_sc_Tro]',noise);

%---------------- SVALBARD -----------------------------------------------
% Initialization
reci_sc_Sva = zeros(3,length(t_vis_Sva));
veci_sc_Sva = zeros(3,length(t_vis_Sva));

for i = 1:length(t_vis_Sva)
    % SGP4 propagation
    tsince = (t_vis_Sva(i) - sat_epoch_et)/60.0; % minutes from TLE epoch
    [~,rteme_sc,vteme_sc] = sgp4(satrec,  tsince);

    % Precession
    ttt = cspice_unitim(t_vis_Sva(i), 'ET', 'TDT')/cspice_jyear()/100;

    % TEME to ECI conversion
    [reci_sc_Sva(:,i), veci_sc_Sva(:,i)] = ...
        teme2eci(rteme_sc, vteme_sc, [0.0;0.0;0.0],  ttt, ddpsi, ddeps);
end

stationName = SVALBARD.name;
noise = SVALBARD.noise;

% SGP4 without noise
[az_Sva, elev_Sva, ~, ~] = antennaPointing(stationName, t_vis_Sva,[reci_sc_Sva; veci_sc_Sva]');
% SGP4 with noise
[az_Sva_meas, elev_Sva_meas, range_Sva_meas, rangeRate_Sva_meas] = antennaPointing(stationName, t_vis_Sva,[reci_sc_Sva; veci_sc_Sva]',noise);
%-------------------------------------------------------------------------

% Plot of the different cases

% New Visibility index SGP4 with errors
i_vis_Kou_meas = elev_Kou_meas > deg2rad(KOUROU.min_elevation);
i_vis_Tro_meas = elev_Tro_meas > deg2rad(TROLL.min_elevation);
i_vis_Sva_meas = elev_Sva_meas > deg2rad(SVALBARD.min_elevation);

% Visibility index with SPG4
i_vis_Kou_sgp4 = elev_Kou > deg2rad(KOUROU.min_elevation);
i_vis_Tro_sgp4 = elev_Tro > deg2rad(TROLL.min_elevation);
i_vis_Sva_sgp4 = elev_Sva > deg2rad(SVALBARD.min_elevation);

figure;
hold on
% % % % Plot SGP4
% plot(az_Kou(i_vis_Kou_sgp4)*cspice_dpr(), elev_Kou(i_vis_Kou_sgp4)*cspice_dpr(),'x','MarkerSize',8,'LineWidth',1.2)
% plot(az_Tro(i_vis_Tro_sgp4)*cspice_dpr(), elev_Tro(i_vis_Tro_sgp4)*cspice_dpr(),'*','MarkerSize',8,'LineWidth',1.2)
% plot(az_Sva(i_vis_Sva_sgp4)*cspice_dpr(), elev_Sva(i_vis_Sva_sgp4)*cspice_dpr(),'+','MarkerSize',8,'LineWidth',1.2)
% hold on

%Plot with error
plot(az_Kou_meas(i_vis_Kou_meas)*cspice_dpr(), elev_Kou_meas(i_vis_Kou_meas)*cspice_dpr(),'x','DisplayName', 'SMOS','MarkerSize',10,'LineWidth',1.2)
plot(az_Tro_meas(i_vis_Tro_meas)*cspice_dpr(), elev_Tro_meas(i_vis_Tro_meas)*cspice_dpr(),'*','DisplayName', 'SMOS','MarkerSize',10,'LineWidth',1.2)
plot(az_Sva_meas(i_vis_Sva_meas)*cspice_dpr(), elev_Sva_meas(i_vis_Sva_meas)*cspice_dpr(),'+','DisplayName', 'SMOS','MarkerSize',10,'LineWidth',1.2)
%Plot Keplerian
plot(azimuth_Kou(i_vis_Kou)*cspice_dpr(), elevation_Kou(i_vis_Kou)*cspice_dpr(),'x','DisplayName', 'SMOS','MarkerSize',10,'LineWidth',1.2)
plot(azimuth_Tro(i_vis_Tro)*cspice_dpr(), elevation_Tro(i_vis_Tro)*cspice_dpr(),'*','DisplayName', 'SMOS','MarkerSize',10,'LineWidth',1.2)
plot(azimuth_Sva(i_vis_Sva)*cspice_dpr(), elevation_Sva(i_vis_Sva)*cspice_dpr(),'+','DisplayName', 'SMOS','MarkerSize',10,'LineWidth',1.2)
grid on;
axis([-180,180,0, 90])
xlabel('Azimuth [deg]')
ylabel('Elevation [deg]')
legend({'Kourou measured values', 'Troll measured values', 'Svalbard measured values', 'Kourou with Keplerian motion', 'Troll with Keplerian motion ','Svalbard with Keplerian motion '}, ...
    'NumColumns', 2);

title('Visibility Windows')


%% --------------2.3 Solve the navigation problem ----------------------

% Initial state
x0 = yy_0(end,:); % initial guess from propagation with tbp from the reference epoch

% Propagation time
t_span = tspan_1;

%----------------------------- Case a ----------------------------------

% Get Station
Stations_a = struct('stationName', {}, 'data', {}, 'visibilityWindow', {}, 'W_m', {}, 'Measurments',{});
% --------- KOUROU ------------
Stations_a(1).stationName = KOUROU.name;
Stations_a(1).data = KOUROU;
Stations_a(1).visibilityWindow = t_vis_Kou;
Stations_a(1).W_m = inv(diag([KOUROU.noise.sigma_az_el /cspice_dpr();KOUROU.noise.sigma_az_el /cspice_dpr() ;KOUROU.noise.sigma_range]));
Stations_a(1).Measurments = [az_Kou_meas, elev_Kou_meas, range_Kou_meas];

% Ode function
ODEfun =  @(t,x) ode_2bp(t,x,satrec.mu);
% Cost function
fun = @(x) costFun(x, t_span, Stations_a,  ODEfun);
% Call lsqnonlin
options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter','StepTolerance',1e-10);

[x_a,resnorm_a,residual_a,exitflag_a,~,~,jacobian_a] = lsqnonlin(fun, x0, [], [], options);

% Covariance computation
Jac_a = full(jacobian_a);
P_ls_a = resnorm_a/(length(residual_a)-length(x0)).*inv(Jac_a.'*Jac_a);

% Linear Mapping
[std_a, std_i] = linearMapping(x_a, P_ls_a, satrec.mu);

% Navigation solution
NavSolution_a.State = x_a;
NavSolution_a.Cov = P_ls_a;
NavSolution_a.Jac = jacobian_a;
NavSolution_a.sigmaPos = sqrt(trace(P_ls_a(1:3,1:3)));
NavSolution_a.sigmaVel = sqrt(trace(P_ls_a(4:6,4:6)));
NavSolution_a.sigmaSMA = std_a;
NavSolution_a.sigmaInc = std_i;


%-------------------------CASE b ----------------------------------------

Stations = struct('stationName', {}, 'data', {}, 'visibilityWindow', {}, 'W_m', {}, 'Measurments',{});

% --------- KOUROU ------------
Stations(1).stationName = KOUROU.name;
Stations(1).data = KOUROU;
Stations(1).visibilityWindow = t_vis_Kou;
Stations(1).W_m = inv(diag([KOUROU.noise.sigma_az_el /cspice_dpr();KOUROU.noise.sigma_az_el /cspice_dpr() ;KOUROU.noise.sigma_range]));
Stations(1).Measurments = [az_Kou_meas, elev_Kou_meas, range_Kou_meas];

% ---------- TROLL -------------
Stations(2).stationName = TROLL.name;
Stations(2).data = TROLL;
Stations(2).visibilityWindow = t_vis_Tro;
Stations(2).W_m = inv(diag([TROLL.noise.sigma_az_el /cspice_dpr();TROLL.noise.sigma_az_el /cspice_dpr() ;TROLL.noise.sigma_range]));
Stations(2).Measurments = [az_Tro_meas, elev_Tro_meas, range_Tro_meas];

% ---------- SVALBARD -------------
Stations(3).stationName = SVALBARD.name;
Stations(3).data = SVALBARD;
Stations(3).visibilityWindow = t_vis_Sva;
Stations(3).W_m = inv(diag([SVALBARD.noise.sigma_az_el /cspice_dpr();SVALBARD.noise.sigma_az_el /cspice_dpr() ;SVALBARD.noise.sigma_range]));
Stations(3).Measurments = [az_Sva_meas, elev_Sva_meas, range_Sva_meas];

% Batch filter with lsqnonlin
ODEfun =  @(t,x) ode_2bp(t,x,satrec.mu);
fun =@(x) costFun(x, t_span, Stations,  ODEfun);

options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter','StepTolerance',1e-10);
[x_b,resnorm_b,residual_b,exitflag_b,~,~,jacobian_b] = lsqnonlin(fun, x0, [], [], options);

% Post-processing
Jac_b = full(jacobian_b);
P_ls_b = resnorm_b/(length(residual_b)-length(x0)).*inv(Jac_b.'*Jac_b);

% Linear Mapping
[std_a, std_i] = linearMapping(x_b, P_ls_b, satrec.mu);

% Navigation solution
NavSolution_b.State = x_b;
NavSolution_b.Cov = P_ls_b;
NavSolution_b.Jac = jacobian_b;
NavSolution_b.sigmaPos = sqrt(trace(P_ls_b(1:3,1:3)));
NavSolution_b.sigmaVel = sqrt(trace(P_ls_b(4:6,4:6)));
NavSolution_b.sigmaSMA = std_a;
NavSolution_b.sigmaInc = std_i;

%------------------------------ Case c ------------------------------------
% J2 
J2 = 0.0010826269;
E_RADII = cspice_bodvrd('Earth','RADII',3);

% Batch filter with lsqnonlin
ODEfun =  @(t,x) ode_2bp_perturbed( t, x, satrec.mu, J2, E_RADII(1));
fun =@(x) costFun(x, t_span, Stations,  ODEfun);

options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter','StepTolerance',1e-10);

[x_c,resnorm_c,residual_c,exitflag_c,~,~,jacobian_c] = lsqnonlin(fun, x0, [], [], options);

% Post-procesing
Jac_c = full(jacobian_c);
P_ls_c = resnorm_c/(length(residual_c)-length(x0)).*inv(Jac_c.'*Jac_c);

% Linear Mapping
[std_a, std_i] = linearMapping(x_c, P_ls_c, satrec.mu);

% Navigation Solution
NavSolution_c.State = x_c;
NavSolution_c.Cov = P_ls_c;
NavSolution_c.Jac = jacobian_c;
NavSolution_c.sigmaPos = sqrt(trace(P_ls_c(1:3,1:3)));
NavSolution_c.sigmaVel = sqrt(trace(P_ls_c(4:6,4:6)));
NavSolution_c.sigmaSMA = std_a;
NavSolution_c.sigmaInc = std_i;




%% ------------------------- 2.4 Trade-off analysis -----------------------

%-------------- CASE 1 KOUROU - TROLL ----------------------------------
% Batch filter with lsqnonlin
ODEfun =  @(t,x) ode_2bp_perturbed( t, x, satrec.mu, J2, E_RADII(1));
fun = @(x) costFun(x, t_span, Stations(1:2),  ODEfun);

options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter','StepTolerance',1e-10);
[x_KT,resnorm_KT,residual_KT,exitflag_KT,~,~,jacobian_KT] = lsqnonlin(fun, x0, [], [], options);

% Post-procesing
Jac_KT = full(jacobian_KT);
P_ls_KT = resnorm_KT/(length(residual_KT)-length(x0)).*inv(Jac_KT.'*Jac_KT);

% Linear Mapping
[std_a, std_i] = linearMapping(x_KT, P_ls_KT, satrec.mu);

% Navigation Solution
NavSolution_KT.State = x_KT;
NavSolution_KT.Cov = P_ls_KT;
NavSolution_KT.Jac = jacobian_KT;
NavSolution_KT.sigmaPos = sqrt(trace(P_ls_KT(1:3,1:3)));
NavSolution_KT.sigmaVel = sqrt(trace(P_ls_KT(4:6,4:6)));
NavSolution_KT.sigmaSMA = std_a;
NavSolution_KT.sigmaInc = std_i;


% --------------   CASE 2 KOUROU - SVALBARD -------------------------
% Batch filter with lsqnonlin
ODEfun =  @(t,x) ode_2bp_perturbed( t, x, satrec.mu, J2, E_RADII(1));
fun = @(x) costFun(x, t_span, Stations(1:2:3),  ODEfun);

options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter','StepTolerance',1e-10);
[x_KS,resnorm_KS,residual_KS,exitflag_KS,~,~,jacobian_KS] = lsqnonlin(fun, x0, [], [], options);

% Post-processing
Jac_KS = full(jacobian_KS);
P_ls_KS = resnorm_KS/(length(residual_KS)-length(x0)).*inv(Jac_KS.'*Jac_KS);

% Linear Mapping
[std_a, std_i] = linearMapping(x_KS, P_ls_KS, satrec.mu);

% Navigation Solution
NavSolution_KS.State = x_KS;
NavSolution_KS.Cov = P_ls_KS;
NavSolution_KS.Jac = jacobian_KS;
NavSolution_KS.sigmaPos = sqrt(trace(P_ls_KS(1:3,1:3)));
NavSolution_KS.sigmaVel = sqrt(trace(P_ls_KS(4:6,4:6)));
NavSolution_KS.sigmaSMA = std_a;
NavSolution_KS.sigmaInc = std_i;


% ------------ CASE 3 TROLL - SVALBARD ------------------------
% Batch filter with lsqnonlin
ODEfun =  @(t,x) ode_2bp_perturbed( t, x, satrec.mu, J2, E_RADII(1));
fun = @(x) costFun(x, t_span, Stations(2:3),  ODEfun);

options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'iter','StepTolerance',1e-10);
[x_TS,resnorm_TS,residual_TS,exitflag_TS,~,~,jacobian_TS] = lsqnonlin(fun, x0, [], [], options);

% Post-processing
Jac_TS = full(jacobian_TS);
P_ls_TS = resnorm_TS/(length(residual_TS)-length(x0)).*inv(Jac_TS.'*Jac_TS);

% Linear Mapping
[std_a, std_i] = linearMapping(x_TS, P_ls_TS, satrec.mu);

% Navigation Solution
NavSolution_TS.State = x_TS;
NavSolution_TS.Cov = P_ls_TS;
NavSolution_TS.Jac = jacobian_TS;
NavSolution_TS.sigmaPos = sqrt(trace(P_ls_TS(1:3,1:3)));
NavSolution_TS.sigmaVel = sqrt(trace(P_ls_TS(4:6,4:6)));
NavSolution_TS.sigmaSMA = std_a;
NavSolution_TS.sigmaInc = std_i;

% Results
disp('KOUROU TROLL')
disp ('Inc [deg] , SMA [km]')
disp(NavSolution_KT.sigmaInc)
disp(NavSolution_KT.sigmaSMA)

disp('KOUROU SVALBARD')
disp ('Inc [deg] , SMA [km]')
disp(NavSolution_KS.sigmaInc)
disp(NavSolution_KS.sigmaSMA)

disp('TROLL SVALBARD')
disp ('Inc [deg] , SMA [km]')
disp(NavSolution_TS.sigmaInc)
disp(NavSolution_TS.sigmaSMA)

%% ---------------------- 2.5 Long Term Analysis --------------------------

% Initialization
y_ref = [reci;veci];
delta_t = 12*3600;  % 12 hours
% forward and backward propagation around the reference tle
et_back = et_ref - delta_t;
et_forw = et_ref + delta_t;
npoints = round((et_ref-et_back)/30.0) + 1;
tspan_back = linspace(et_ref, et_back, npoints);
tspan_forw = linspace(et_ref, et_forw, npoints);
% Options
options = odeset('RelTol',1e-13,'AbsTol',1e-13);

% Propagation from starting time to final time
[tt_forw,xx_forw] = ode113(@(t,y) ode_2bp(t,y,satrec.mu), tspan_forw, y_ref,options);

% Propagation from reference to starting time
[tt_back,xx_back] = ode113(@(t,y) ode_2bp(t,y,satrec.mu), tspan_back, y_ref,options);

% Assembling
tt = [flip(tt_back);tt_forw(2:end)];
yy_ECI = [flip(xx_back,1);xx_forw(2:end,:)];

% Visibility
stationName = KOUROU.name;
[azimuth_Kou, elevation_Kou, range_Kou, rangeRate_Kou] = antennaPointing(stationName, tt(1:end), yy_ECI(1:end,:));

stationName = TROLL.name;
[azimuth_Tro, elevation_Tro, range_Tro, rangeRate_Tro] = antennaPointing(stationName, tt, yy_ECI);

stationName = SVALBARD.name;
[azimuth_Sva, elevation_Sva, range_Sva, rangeRate_Sva] = antennaPointing(stationName, tt(1:end), yy_ECI(1:end,:));

% Visibility indices
i_vis_Kou = elevation_Kou > deg2rad(KOUROU.min_elevation);
i_vis_Tro = elevation_Tro > deg2rad(TROLL.min_elevation);
i_vis_Sva = elevation_Sva > deg2rad(SVALBARD.min_elevation);

figure
plot(azimuth_Kou(i_vis_Kou)*cspice_dpr(), elevation_Kou(i_vis_Kou)*cspice_dpr(),'*','DisplayName', 'SMOS','MarkerSize',8,'LineWidth',1.2)
hold on
plot(azimuth_Tro(i_vis_Tro)*cspice_dpr(), elevation_Tro(i_vis_Tro)*cspice_dpr(),'*','DisplayName', 'SMOS','MarkerSize',8,'LineWidth',1.2)
plot(azimuth_Sva(i_vis_Sva)*cspice_dpr(), elevation_Sva(i_vis_Sva)*cspice_dpr(),'*','DisplayName', 'SMOS','MarkerSize',8,'LineWidth',1.2)
grid on;
axis([-180,180,0, 90])
xlabel('Azimuth [deg]')
ylabel('Elevation [deg]')
legend('Kourou (60 s)','Troll (30 s)','Svalbard (60 s)')
title('Visibility Windows')

% Select only visible interval
elevation_Kou_filt = elevation_Kou; 
elevation_Kou_filt(elevation_Kou < deg2rad(KOUROU.min_elevation)) = NaN;

elevation_Tro_filt = elevation_Tro; 
elevation_Tro_filt(elevation_Tro < deg2rad(TROLL.min_elevation)) = NaN;

elevation_Sva_filt = elevation_Sva; 
elevation_Sva_filt(elevation_Sva < deg2rad(SVALBARD.min_elevation)) = NaN;

%
figure
plot(tt/cspice_spd,elevation_Kou_filt*cspice_dpr(),'-','DisplayName', 'SMOS','MarkerSize',8,'LineWidth',1.2)
hold on
plot(tt/cspice_spd,elevation_Tro_filt*cspice_dpr(),'-','DisplayName', 'SMOS','MarkerSize',8,'LineWidth',1.2)
plot(tt/cspice_spd,elevation_Sva_filt*cspice_dpr(),'-','DisplayName', 'SMOS','MarkerSize',8,'LineWidth',1.2)
grid on;
ylim([0, 90])
xlabel('Epoch [MJD2000]')
ylabel('Elevation [deg]')
legend('Kourou ','Troll ','Svalbard ')
title('Vibility Windows')
tick_values = linspace(tt(1)/cspice_spd, tt(end)/cspice_spd, 4); % 4 valori equidistanti

% Imposta i tick sull'asse X
xticks(tick_values);


cspice_kclear()

%% FUNCTIONS

function dy = ode_2bp( ~, y, mu )
%----------------------------------------------------------------
% ode_2bp ODE system for the two-body problem (Keplerian motion)
%
% INPUT:
% t         [1] Time (can be omitted, as the system is autonomous) 
% y         [6x1] State of the body ( rx, ry, rz, vx, vy, vz ) 
% mu        [1] Gravitational parameter of the primary 
%
% OUTPUT:
% dy        [6x1] Derivative of the state 
% -------------------------------------------------------------------------

% Position and velocity
r = y(1:3);
v = y(4:6);
% Distance from the primary
rnorm = norm(r);
% Set the derivatives of the state
dy = [ v ; (-mu/rnorm^3)*r ];
end

function dy = ode_2bp_perturbed( t, y,mu_E, J2, R_E)
%----------------------------------------------------------------
% ode_2bp ODE system for the two-body problem with the
% J2 perturbation at each time
%
% INPUT:
% t         [1] Time (can be omitted, as the system is autonomous) 
% y         [6x1] State of the body ( rx, ry, rz, vx, vy, vz ) 
% mu_E      [1] Gravitational parameter of the primary 
% J2        [1] Earth's second spherical armonic coefficent
% R_E       [1] Earrth's equatorial radius
%
% OUTPUT:
% dy        [6x1] Derivative of the state 
% -------------------------------------------------------------------------


rECI = y(1:3);
vECI = y(4:6);
% Distance from the primary
r = norm(rECI);

% Rotation to ECEF
rotM= cspice_pxform('J2000', 'ITRF93', t);
rECEF = rotM * rECI ;
rECEFnorm   = norm(rECEF) ;

% J2 effect
J2_coeff = (3/2 * J2 * mu_E * R_E^2/rECEFnorm^5);

aJ2ECEF =   [J2_coeff*rECEF(1)*(5*rECEF(3)^2/rECEFnorm^2 - 1); 
             J2_coeff*rECEF(2)*(5*rECEF(3)^2/rECEFnorm^2 - 1);
             J2_coeff*rECEF(3)*(5*rECEF(3)^2/rECEFnorm^2 - 3)];
aJ2ECI = rotM' * aJ2ECEF ;

dy = [vECI ; (-mu_E/r^3)*rECI + aJ2ECI];
end

function [azimuth, elevation, range, range_rate] = antennaPointing(stationName, et_vec, yy_ECI,noise)
%----------------------------------------------------------------
% Evaluate the measurments of azimuth elevation and range from a given
% ground reference with available kernel from a second body
% set of cordinates over a time window
%
% INPUT:
% stationName  [str] Name of the ground reference
% yy_ECI       [nx6] State of the ( rx, ry, rz, vx, vy, vz ) [km, km/s]
% et_vec       [nx1] Time interval 
% noise        [3x3] Diagonal noise matrix
%
% OUTPUT:
% azimuth      [nx1] Azimuth vector [rad]
% elevation    [nx1] Elevation vector [rad]
% range        [nx1] Range vector [km]
% range_rate   [nx1] Range rate vector [km/s]
% -------------------------------------------------------------------------

% Topocentric rf for the ground station
station_TOPO = [stationName, '_TOPO'];

% Initialization
azimuth = zeros(length(et_vec),1);
elevation = zeros(length(et_vec),1);
range = zeros(length(et_vec),1);
range_rate = zeros(length(et_vec),1);

for j = 1:length(et_vec)
    % Compute station intial position in ECI
    rv_station_eci = cspice_spkezr(stationName, et_vec(j), 'J2000', 'NONE', 'EARTH');
    
    % Compute station-satellite vector in ECI
    rv_station_sat_eci = yy_ECI(j,:)' - rv_station_eci;
    
    % Transformation from ECI to topocentric frame
    ROT_ECI2TOPO = cspice_sxform('J2000', station_TOPO, et_vec(j));

    % Convert state into topocentric frame
    rv_station_sat_topo = ROT_ECI2TOPO*rv_station_sat_eci;

    % Compute range, azimuth and elevation using cspice_xfmsta
    rll_station_sat = cspice_xfmsta(rv_station_sat_topo,'RECTANGULAR','LATITUDINAL','EARTH');
    
    range(j)   = rll_station_sat(1);   % [km]
    range_rate(j) = dot(rv_station_sat_topo(4:6),rv_station_sat_topo(1:3)) / range(j); % [km/s]
    azimuth(j) = rll_station_sat(2);   % [rad]
    elevation(j) = rll_station_sat(3); % [rad]
   
end

    if nargin == 4

       % Standard Deviations
       sigma_az_el = deg2rad(noise.sigma_az_el);
       sigma_range = noise.sigma_range;
       
       % Gaussian Errors
       range_old = range;
       range = mvnrnd(range,sigma_range^2);
       range_rate = range_rate .* range_old ./ range; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       azimuth = mvnrnd(azimuth,sigma_az_el^2);
       elevation = mvnrnd(elevation,sigma_az_el^2);

    end

end

function residual = costFun(x, t_span, Stations,  ODEfun)
%----------------------------------------------------------------
% Cost function for the lsqnonlin matlab function in order to implement a
% batch filter. return the residuals beyween the estimated measurments and
% a set of real measurments from different ground stations.
%
% INPUT:
% t_span    [nx1] Time vector 
% x         [6x1] Initial state of the body
% ODEfun    [func] Dynamics model for the integration
%
% OUTPUT:
% residual  [mxl] Weighted residual matrix
% -------------------------------------------------------------------------
% Get number of stations
nStations = length(Stations);

% Propagation
options = odeset('Reltol', 1.e-13, 'Abstol', 1.e-13);
[~, xx_ECI] = ode113(ODEfun, t_span, x, options);

% Pre-allocation of the residual
n_meas_tot = sum(arrayfun(@(s) size(s.Measurments, 1), Stations));
residual = zeros(n_meas_tot, size(Stations(1).Measurments, 2));

%
i_old = 1;
    for i = 1:nStations
           
        % Station i data
        t_vis = Stations(i).visibilityWindow;
        W_m = Stations(i).W_m;
        meas_real = Stations(i).Measurments;
    
        % Get Measurment time position
        [isMeasurTime, ~] = ismember(t_span, t_vis);
        meas_indices = find(isMeasurTime);
        % Predicted measurements
        [azimuth, elevation, range, ~] = antennaPointing(Stations(i).stationName, t_span(meas_indices), xx_ECI(meas_indices,:));
        meas_pred = [azimuth, elevation, range];
    
        % Residuals station i
        diff_meas_ang = (W_m(1:2,1:2) * angdiff(meas_pred(:,1:2)', meas_real(:,1:2)'))';
        diff_meas_lon = W_m(3,3) * (meas_pred(:,3) - meas_real(:,3));
        res = [diff_meas_ang, diff_meas_lon];
    
        % Total residual
        i_new = i_old + size(meas_real, 1) - 1;
        residual(i_old:i_new, :) = res;
        i_old = i_new + 1;
    end
end

function [std_a, std_i] = linearMapping(stateVec, P, mu)
% -----------------------------------------------------------------
% Compute the linear mapping from cartesian state to keplerian elements, to retrieve the standard deviation of semimajor
% axis and inclination.
%
% INPUT:
%   stateVec   [6x1] Cartesian state vector [x, y, z, vx, vy, vz]
%   P          [6x6] Covariance matrix of the cartesian state vector
%   mu         [1x1] Gravitational parameter of the central body [km^3/s^2].
%
% OUTPUT:
%   std_a      [1x1] standard deviation of semimajor axis [km]
%   std_i      [1x1] standard deviation of inclination [deg]
% ----------------------------------------------------------------
    
% Initialization
Jac = zeros(6, 6);
rr = stateVec(1:3)';
vv = stateVec(4:6)';
 
% position
for k = 1:3
    eps_r = zeros(3,1);
    eps_r(k) = sqrt(eps) * max(1, abs(rr(k)));

    % Forward perturbation
    [a_fin, e_fin, i_fin, OM_fin, om_fin, th_fin] = car2kep(rr + eps_r, vv, mu);
    [a_in, e_in, i_in, OM_in, om_in, th_in] = car2kep(rr, vv, mu);
    % Jacobian from Cartesian to Keplerian elements.
    % Numerical derivative with finite difference
    Jac(:, k) = ([a_fin, e_fin, i_fin, OM_fin, om_fin, th_fin]' - [a_in, e_in, i_in, OM_in, om_in, th_in]') ./ eps_r(k);
end

% velocity
for k = 1:3
    eps_v = zeros(3,1);
    eps_v(k) = sqrt(eps) * max(1, abs(vv(k)));

    % Forward and backward perturbation
    [a_fin, e_fin, i_fin, OMf, om_fin, th_fin] = car2kep(rr, vv + eps_v, mu);
    [a_in, e_in, i_in, OM_in, om_in, th_in] = car2kep(rr, vv, mu);
    % Jacobian from Cartesian to Keplerian elements
    % Numerical derivative with finite difference
    Jac(:, k+3) = ([a_fin, e_fin, i_fin, OMf, om_fin, th_fin]'- [a_in, e_in, i_in, OM_in, om_in, th_in]') ./ eps_v(k);
end

% Compute the covariance
P_kep = Jac * P * Jac';

std_a = sqrt(P_kep(1,1)); % Standard deviation of semimajor axis [km]
std_i = rad2deg(sqrt(P_kep(3,3))); % Standard deviation of inclination [deg]
end

function [a, e, i, OM, om, th,N] = car2kep(r, v, mu)
% car2par_rad Converts Cartesian coordinates to Keplerian parameters (Output in radians)
%
% INPUT:
%   r       [3x1] Position vector [km]
%   v       [3x1] Velocity vector [km/s]
%   mu      [1x1] Gravitational parameter [km^3/s^2]
%
% OUTPUT:
%   a       [1x1] Semi-major axis [km]
%   e       [1x1] Eccentricity [-]
%   i       [1x1] Inclination [rad]
%   OM      [1x1] Right Ascension of the Ascending Node (RAAN) [rad]
%   om      [1x1] Argument of periapsis [rad]
%   th      [1x1] True anomaly [rad]

I = [1 0 0]';
J = [0 1 0]';
K = [0 0 1]';

r_norm = norm(r);
v_norm = norm(v);

h = cross(r,v);
h_norm = norm(h);


i = acos(dot(h,K)/h_norm);


e = 1/mu * ((v_norm^2 - mu/r_norm)*r - (dot(r,v))*v);
e_norm = norm(e);


Eps = 1/2 * v_norm^2 - mu/r_norm;
a = -mu/(2*Eps);

N = cross(K,h);
N_norm = norm(N);


if dot(N,J) >= 0
 OM = acos(dot(N,I)/N_norm);
else
    OM = 2*pi - acos(dot(N,I)/N_norm);
end
if i == 0
    OM = 0;
    N = [1;0;0];
    N_norm = norm(N);
end

if dot(e,K) >= 0
    om = acos(dot(N,e)/(N_norm * e_norm));
else
    om = 2*pi -acos(dot(N,e)/(N_norm * e_norm));
end

if dot(v,r)/r_norm >=0
    th = acos(dot(r,e)/(r_norm * e_norm));
else
    th = 2*pi - acos(dot(r,e)/(r_norm * e_norm));
end

e = e_norm;
end

function plot_set
%-------------------------------------
% Plots settings
%--------------------------------------
% Interpreter:
set(0, 'defaultTextInterpreter', 'latex')
set(0, 'defaultLegendInterpreter', 'latex')
set(0, 'defaultAxesTickLabelInterpreter', 'latex')
% % Setting Legends:
% set(0, 'defaultLegendLocation','southwest');
% set(0, 'defaultLegendOrientation', 'vertical');
% set(0, 'defaultLegendFontSize', 12);
% Setting Axes:
set(0, 'defaultAxesXMinorGrid', 'on');
set(0,'defaultAxesYMinorGrid','on');
set(0, 'defaultAxesFontSize', 15);

end