% Spacecraft Guidance and Navigation (2024/2025)
% Assignment # 2, Exercise 3
% Author: Davide Bellini

clearvars; close all; clc;
cspice_kclear()


plot_set(); % Plot settings
% rng setted to default to provides always the same results, if necessary
% remove it
rng default
% Upload kernels
% Load kernels
cspice_furnsh('assignment02.tm');


%% DATA
% Initial conditions
r0 = [4307.844185282820, -1317.980749248651, 2109.210101634011]'; % [km]
v0 = [-0.110997301537882, -0.509392750828585, 0.815198807994189]'; %[km/s]
% Initial and Final time  [UTC] 
initTime = '2024-11-18T16:30:00.000'; 
finalTime = '2024-11-18T20:30:00.000';
% Measurements noise 
sigma_p = 100 * 1e-3; % [km]
% Covariance P0 [km2, km2/s2, rad2] 
P0 = diag([10,1,1,0.001,0.001,0.001,0.00001,0.00001]);

% Lander
Moonlander.LanderName = 'MOONLANDER';
Moonlander.coord.lat = 78;       % [deg]
Moonlander.coord.lon = 15;       % [deg]
Moonlander.coord.alt = 0;        % [km]
Moonlander.minElevation = 0;     % [deg]

% UT parameter
alpha = 0.01;
beta = 2;


%% ----------------------- 3.1 Check the visibility window ----------------

% Get Initial and Final et 
et_0 = cspice_str2et(initTime);
et_f = cspice_str2et(finalTime);

% Time span each 30 sec
dt = 30; % [s]
tspan = et_0 : dt : et_f;

% Get moon gravity parameter
mu_M = cspice_bodvrd('Moon','GM',1);   % mu Moon

% Propagation
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
[tt,xx_MCI] = ode113(@(t,y) ode_2bp(t,y,mu_M), tspan, [r0;v0],options);

% Radians conversion
lat = Moonlander.coord.lat * cspice_rpd();   
lon = Moonlander.coord.lon * cspice_rpd();  
alt = Moonlander.coord.alt;

% Moon Radii
Mradii = cspice_bodvrd('MOON', 'RADII', 3); 
Mradii_eq= Mradii(1);              
vLand_IAU = [0;0;0];
xxLand_MCI = zeros(length(tspan),6);


for i = 1:length(tspan)
    
    et = tspan(i);
    % Conversion (MCMF)
    rLand_IAU = cspice_latrec(Mradii_eq + alt, lon, lat);
    % Rotation from MCMF to J2000 
    rotMat = cspice_sxform('IAU_MOON', 'J2000', et);

    xxLand_MCI(i,:) = rotMat* [rLand_IAU;vLand_IAU];
end

% Plot orbit %%%%%%%%%%%%%%%%%%%%%%%  togliere immagine %%%%%%%%%%%%%%%%
[X, Y, Z] = sphere(50);
%texture = imread('Moon.jpg');
figure()
hold on
grid on
plot3(xx_MCI(:,1), xx_MCI(:,2), xx_MCI(:,3),'LineWidth',2);
plot3(xx_MCI(1,1), xx_MCI(1,2), xx_MCI(1,3), 'Marker','d','MarkerSize',10,'LineWidth',1.5,'LineStyle','none');
plot3(xx_MCI(end,1), xx_MCI(end,2), xx_MCI(end,3), 'Marker','d','MarkerSize',10,'LineWidth',1.5,'LineStyle','none');
plot3(xxLand_MCI(1,1)+5, xxLand_MCI(1,2)+5, xxLand_MCI(1,3)+5, 'Marker','o','MarkerSize',6,...
       'LineWidth',1,'LineStyle','none','Color','k','MarkerFaceColor',[0.4, 0.9, 0.4]);
hSurface = surf(X*Mradii_eq, Y*Mradii_eq, Z*Mradii_eq,'EdgeColor','none');
colormap('bone'); 
%set(gca, 'CLim', [0, 255]); 
%surface = findobj(gca, 'Type', 'Surface'); 
%set(hSurface, 'FaceColor', 'texturemap', 'CData', texture); 
legend('Trajectory', 'Initial Position', 'Final Position','Lander Initial position', 'Moon')
axis equal;
xlabel('Position X [km]');
ylabel('Position Y [km]');
zlabel('Position Z [km]');
%------------------------------------------------------------------------

% Visibility window
[~, elevation, ~] = visibility(xxLand_MCI, xx_MCI, tspan);

elevation_filt = elevation; 
elevation_filt(elevation < deg2rad(Moonlander.minElevation)) = NaN;

% Plot visibility 
figure
plot(tspan/cspice_spd, elevation_filt*cspice_dpr(),'-','LineWidth',1.4,'Color',[0, 0.4470, 0.7410])
hold on
grid on
plot(tspan/cspice_spd,Moonlander.minElevation*ones(length(elevation)),'--','LineWidth',1.4,'Color',[0.8500, 0.3250, 0.0980])
ylim([-90 90])
xlim([(tspan(1)-100)/cspice_spd, (tspan(end)+100)/cspice_spd])
xlabel('Epoch [MJD2000]')
ylabel('elevation [deg]')
legend('Satellite elevation','Minimum elevation for the visibility','Location','southeast')
title('Visibility window')

% Chech visibility
if any(isnan(elevation_filt))
    disp('Relative visibility not always guarntee')
else
    disp('Relative visibility for the entire time interval ')
end

%% ---------------------- 3.2 Simulate measurements-----------------------

% Error on the satellite position vector
rrMCI_err = mvnrnd(xx_MCI(:,1:3),diag([sigma_p^2,sigma_p^2,sigma_p^2]));
% Full state 
xxMCI_err = [rrMCI_err,xx_MCI(:,4:6)];

% Range measurment
% Without Noise
[azimuth, elevation, range] = measurments(Moonlander.LanderName, tspan, xx_MCI);
% With Noise
[azimuth_meas, elevation_meas, range_meas] = measurments(Moonlander.LanderName, tspan, xx_MCI, sigma_p);
% The range measurment is affected only by the noise on the range
% measurment, the range measurment from the lunar navigation center will be
% affected by both error but they are referred to 2 different measurment.

% New visibility
elevMeas_filt = elevation_meas; 
elevMeas_filt(elevation_meas < deg2rad(Moonlander.minElevation)) = NaN;

% Plot visibility window with measurment
figure
hold on
plot(tspan/cspice_spd,elevMeas_filt*cspice_dpr(),'MarkerSize',4,'LineWidth',0.8,'Marker','+','LineStyle','none')
plot(tspan/cspice_spd,Moonlander.minElevation*ones(length(elevation_meas)),'--')
plot(tspan/cspice_spd,elevation_filt*cspice_dpr(),'MarkerSize',4,'LineWidth',0.8,'Marker','x','LineStyle','none')
grid on
ylim([-90 90])
xlim([(tspan(1)-100)/cspice_spd, (tspan(end)+100)/cspice_spd])
xlabel('Epoch [MJD2000]')
ylabel('Elevation [deg]')
legend('Satellite elevation measured','Minimum elevation for the visibility','Satellite elevation','Location','southeast')
title('Visibility with Real Measurments')
% Check visibility
if any(isnan(elevMeas_filt))
    disp('Relative visibility not always guarntee')
else
    disp('Relative visibility for the entire time interval ')
end

colErr = [0.8500, 0.3250, 0.0980];  % Color for the error
colBound = [0, 0.4470, 0.7410];     % Color for the std
% Graphical solution of the errors
figure

subplot(2,2,1) % absolute error on position x-component
plot(tspan/cspice_spd,xx_MCI(:,1)-xxMCI_err(:,1),'Color',colErr,'LineWidth',1)
hold on
grid on
plot(tspan/cspice_spd,ones(length(tt))*3*sigma_p,'--','Color',colBound,'LineWidth',1)
plot(tspan/cspice_spd,-ones(length(tt))*3*sigma_p,'--','Color',colBound,'LineWidth',1)
ylim([-0.45 0.45])
ylabel('Absolute Error [km]')
xlim([(tspan(1)-50)/cspice_spd, (tspan(end)+50)/cspice_spd])
xlabel('Epoch [MJD2000]')
title('x-coordinate')

subplot(2,2,2) % absolute error on position y-component
plot(tspan/cspice_spd,xx_MCI(:,2)-xxMCI_err(:,2),'Color',colErr,'LineWidth',1)
hold on
grid on
plot(tspan/cspice_spd,ones(length(tt))*3*sigma_p,'--','Color',colBound,'LineWidth',1)
plot(tspan/cspice_spd,-ones(length(tt))*3*sigma_p,'--','Color',colBound,'LineWidth',1)
ylim([-0.45 0.45])
ylabel('Absolute Error [km]')
xlim([(tspan(1)-50)/cspice_spd, (tspan(end)+50)/cspice_spd])
xlabel('Epoch [MJD2000]')
title('y-coordinate')

subplot(2,2,3)
plot(tspan/cspice_spd,xx_MCI(:,3)-xxMCI_err(:,3),'Color',colErr,'LineWidth',1)
hold on
grid on
plot(tspan/cspice_spd,ones(length(tt))*3*sigma_p,'--','Color',colBound,'LineWidth',1)
plot(tspan/cspice_spd,-ones(length(tt))*3*sigma_p,'--','Color',colBound,'LineWidth',1)
ylim([-0.45 0.45])
ylabel('Absolute Error [km]')
xlim([(tspan(1)-50)/cspice_spd, (tspan(end)+50)/cspice_spd])
xlabel('Epoch [MJD2000]')
title('z-coordinate')

subplot(2,2,4)
plot(tspan/cspice_spd,range - range_meas,'Color',colErr,'LineWidth',1)
hold on
grid on
plot(tspan/cspice_spd,ones(length(tt))*3*sigma_p,'--','Color',colBound,'LineWidth',1)
plot(tspan/cspice_spd,-ones(length(tt))*3*sigma_p,'--','Color',colBound,'LineWidth',1)
ylim([-0.45 0.45])
ylabel('Absolute Error [km]')
xlim([(tspan(1)-50)/cspice_spd, (tspan(end)+50)/cspice_spd])
xlabel('Epoch [MJD2000]')
grid on
title('Range')


%% -------------------3.3 Estimate the lunar orbiter absolute state ------------

% Initialization
k = 0;
n = size(xx_MCI,2);
lambda = alpha^2*(n + k) - n; 

% Weights
W_m0 = lambda / (lambda + n);
W_c0 = lambda / (lambda + n) + (1 - alpha^2 + beta);
Wcm = 1/ (2*(n+lambda));

Weight.covariance =  [W_c0; Wcm*ones(2*n,1)]';
Weight.mean = [W_m0; Wcm*ones(2*n,1)]';

% Noise matrix
sigma_mat = diag([sigma_p^2,sigma_p^2,sigma_p^2]);
% Measurments
%yy_meas = xxMCI_err;
measurment = xxMCI_err(:,1:3);
% Initial condition perturbation
xx_0 = mvnrnd(xx_MCI(1,:),P0(1:6,1:6));
 
% UKF
[xx, Px] = UKF(xx_0, P0(1:6,1:6), Weight, lambda, tspan, mu_M, measurment, sigma_mat);


% Postion norm and Standard deviation
% Norm
vvNorm_diff = vecnorm((xx(:,4:6) - xx_MCI(:,4:6))')'; 
rrNorm_diff = vecnorm((xx(:,1:3) - xx_MCI(:,1:3))')';
% Standard deviations position
std_pos = zeros(1, length(tspan)); 
std_vel = zeros(1, length(tspan));
for i = 1:length(tspan)
    std_pos(i) = sum(diag(sqrtm(Px(1:3,1:3, i))));
    std_vel(i) = sum(diag(sqrtm(Px(4:6,4:6, i))));
end

% Plots
figure % Position
semilogy(tspan/cspice_spd,rrNorm_diff,'LineWidth',1.5)
grid on
hold on
semilogy(tspan/cspice_spd,3*std_pos,'LineWidth',1.5)
ylabel('Position Error [km]')
xlim([(tspan(1)-50)/cspice_spd, (tspan(end)+50)/cspice_spd])
xlabel('Epoch [MJD2000]')
legend('Position Error', '3$\sigma$')
title('Position norm error and asssociated Standard deviation')

figure % Velocity
semilogy(tspan/cspice_spd,vvNorm_diff,'LineWidth',1.5)
grid on
hold on
semilogy(tspan/cspice_spd,3*std_vel,'LineWidth',1.5)
ylabel('Velocity  Error[km]')
xlim([(tspan(1)-50)/cspice_spd, (tspan(end)+50)/cspice_spd])
xlabel('Epoch [MJD2000]')
legend('Velocity Error', '3$\sigma$')
title('Velocity error and asssociated Standard deviation')



figure; 
hold on
plot(tspan/cspice_spd,xx_MCI(:,1)-xx(:,1), 'LineWidth', 1)
plot(tspan/cspice_spd,xx_MCI(:,2)-xx(:,2), 'LineWidth', 1)
plot(tspan/cspice_spd,xx_MCI(:,3)-xx(:,3), 'LineWidth', 1)

plot(tspan/cspice_spd,3*std_pos, '--', 'Color', 'k', 'LineWidth', 1)
plot(tspan/cspice_spd,-3*std_pos, '--', 'Color', 'k', 'LineWidth', 1)
ylim([-0.5 0.5])
grid on
ylabel('Position Error [km]')
xlim([(tspan(1)-50)/cspice_spd, (tspan(end)+50)/cspice_spd])
xlabel('Epoch [MJD2000]')
legend('x-coordinate','y-coordinate','z-coordinate','Position 3$\sigma$')
title('Position Error')

%% ---------------------- 3.4 Estimate the lunar lander coordinates -----------------------
% Initialization
state = [xx_MCI(1,:), Moonlander.coord.lon /cspice_dpr(), Moonlander.coord.lat /cspice_dpr()];
k = 0;
n = length(state);
lambda = alpha^2*(n + k) - n;

 % Weights
    W_m0 = lambda / (lambda + n);
    W_c0 = lambda / (lambda + n) + (1 - alpha^2 + beta);
    Wcm = 1/ (2*(n+lambda));

    Weight.covariance =  [W_c0; Wcm*ones(2*n,1)]';
    Weight.mean = [W_m0; Wcm*ones(2*n,1)]';

% Noise matrix
sigma_mat = diag([sigma_p^2,sigma_p^2,sigma_p^2, sigma_p^2]);

% Measurments
measurment = [measurment, range_meas];

% Initial condition
xx_0 = mvnrnd(state,P0);

% UKF
[xx_2, Px_2] = UKF(xx_0, P0, Weight, lambda, tspan, mu_M, measurment, sigma_mat, Mradii_eq);


% Lander coordinates with SPICE
xMoonlanderIau =   cspice_spkezr(Moonlander.LanderName, tspan, 'IAU_MOON', 'NONE', 'MOON');
[~, lon, lat] = cspice_reclat(xMoonlanderIau(1:3,:));

% Norm
rrNorm_diff_2 = vecnorm((xx_2(:,1:3) - xx_MCI(:,1:3))')';
vvNorm_diff_2 = vecnorm((xx_2(:,4:6) - xx_MCI(:,4:6))')';
err_lon = angdiff(lon', xx_2(:,7))*cspice_dpr;
err_lat = angdiff(lat', xx_2(:,8))*cspice_dpr;

% Standard deviations 
std_pos_2 = zeros(1, length(tspan)); 
std_vel_2 = zeros(1, length(tspan)); 
std_lon = zeros(1, length(tspan));
std_lat = zeros(1, length(tspan));

for i = 1:length(tspan)
    % Position and Velocity std
    std_pos_2(i) = sum(diag(sqrtm(Px_2(1:3,1:3, i))));
    std_vel_2(i) = sum(diag(sqrtm(Px_2(4:6,4:6, i))));
    % Latitude and lon
    std_lon(i) = sqrt(Px_2(7,7,i)) * cspice_dpr;
    std_lat(i) = sqrt(Px_2(8,8,i)) * cspice_dpr;
   
end

% Position
figure
semilogy(tspan/cspice_spd,rrNorm_diff,'LineWidth',1.1,'LineStyle','-.')
grid on
hold on
semilogy(tspan/cspice_spd,3*std_pos,'LineWidth',1.1,'LineStyle','-.')
semilogy(tspan/cspice_spd,rrNorm_diff_2,'LineWidth',1.1)
grid on
hold on
semilogy(tspan/cspice_spd,3*std_pos_2,'LineWidth',1.1)
legend('Position Error: only position measurments','3$\sigma$ Position error: only position measurments ',...
       'Position Error: adding range measurment','3$\sigma$ Position error: adding range measurment ')
ylabel('Position Error [km]')
xlim([tspan(1)/cspice_spd, tspan(end)/cspice_spd])
xlabel('Epoch [MJD2000]')

% Velocity
figure
semilogy(tspan/cspice_spd,vvNorm_diff,'LineWidth',1.1,'LineStyle','-.')
grid on
hold on
semilogy(tspan/cspice_spd,3*std_vel,'LineWidth',1.1,'LineStyle','-.')
semilogy(tspan/cspice_spd,vvNorm_diff_2,'LineWidth',1.1)
grid on
hold on
semilogy(tspan/cspice_spd,3*std_vel_2,'LineWidth',1.1)
legend('Velocity Error: only position measurments','3$\sigma$ Velocity error: only position measurments ',...
       'Velocity Error: adding range measurment','3$\sigma$ Velocity error: adding range measurment ')
ylabel('Velocity Error [km]')
xlim([tspan(1)/cspice_spd, tspan(end)/cspice_spd])
xlabel('Epoch [MJD2000]')


% Latitude Error
figure
semilogy(tspan / cspice_spd, abs(err_lat), 'LineWidth', 1.7, 'DisplayName', 'Error'); hold on;
semilogy(tspan/ cspice_spd, 3*abs(std_lat), 'LineWidth', 1.7, 'DisplayName', '3$\sigma$');
xlim([tspan(1)/cspice_spd, tspan(end)/cspice_spd])
xlabel('Epoch [MJD2000]');
ylabel('Latitude Error [deg]');
title('Latitude Error and $3\sigma$');
legend('Location', 'best');
grid on;

%  Longitude Error
figure
semilogy(tspan / cspice_spd, abs(err_lon), 'LineWidth', 1.7, 'DisplayName', ' Error'); hold on;
semilogy(tspan / cspice_spd, 3*abs(std_lon), 'LineWidth', 1.7, 'DisplayName', '3$\sigma$');
xlim([tspan(1)/cspice_spd, tspan(end)/cspice_spd])
xlabel('Epoch [MJD2000]');
ylabel('Longitude Error [deg]');
title('Longitude Error and $3\sigma$');
legend('Location', 'best');
grid on;

%
cspice_kclear()

%% FUNCTIONS

%------------------------------------------------------------------------
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

%------------------------------------------------------------------------
function [azimuth, elevation, range] = visibility(xxLander_MCI, xx_MCI, et_vec)
%----------------------------------------------------------------
% Evaluate the measurments of azimuth elevation and range from two given
% set of cordinates over a time window
%
% INPUT:
% xxLander_MCI [nx6] State of the reference on ground [km,km/s]
% xx_MCI       [nx6] State of the ( rx, ry, rz, vx, vy, vz ) [km, km/s]
% et_vec       [nx1] Time interval 
%
% OUTPUT:
% azimuth      [nx1] Azimuth vector [rad]
% elevation    [nx1] Elevation vector [rad]
% range        [nx1] Range vector [km]
% -------------------------------------------------------------------------

    % Initialization
    azimuth = zeros(length(et_vec),1);
    elevation = zeros(length(et_vec),1);
    range = zeros(length(et_vec),1);
            
    % Compute Lander-satellite vector in MCI
    rv_land_sat_MCI = xx_MCI - xxLander_MCI;
    
    for j = 1:length(et_vec)
        
        % Transformation from ECI to topocentric frame
        ROT_MCI2TOPO = cspice_sxform('J2000', 'MOONLANDER_TOPO', et_vec(j));
    
        % Convert state into topocentric frame
        rv_land_sat_topo = ROT_MCI2TOPO*rv_land_sat_MCI(j,:)';
    
        % Compute range, azimuth and elevation using cspice_xfmsta
        rll_land_sat = cspice_xfmsta(rv_land_sat_topo,'RECTANGULAR','LATITUDINAL','MOON');
        
        range(j)   = rll_land_sat(1);   % [km]
        azimuth(j) = rll_land_sat(2);   % [rad]
        elevation(j) = rll_land_sat(3); % [rad]
   
    end

    

end

%------------------------------------------------------------------------
function [azimuth, elevation, range] = measurments(landerName, et_vec, xx_MCI,sigma)
%----------------------------------------------------------------
% Evaluate the measurments of azimuth elevation and range from a given
% ground reference with available kernel from a second body
% set of cordinates over a time window
%
% INPUT:
% landerName   [str] Name of the ground reference
% xx_MCI       [nx6] State of the ( rx, ry, rz, vx, vy, vz ) [km, km/s]
% et_vec       [nx1] Time interval 
% sigma        [1] Range noise standard deviation
% OUTPUT:
% azimuth      [nx1] Azimuth vector [rad]
% elevation    [nx1] Elevation vector [rad]
% range        [nx1] Range vector [km]
% -------------------------------------------------------------------------

% Topocentric rf for the ground station
lander_TOPO = [landerName, '_TOPO'];

% Initialization
azimuth = zeros(length(et_vec),1);
elevation = zeros(length(et_vec),1);
range = zeros(length(et_vec),1);

for j = 1:length(et_vec)
    % Compute station intial position in MCI
    rv_lander_MCI = cspice_spkezr(landerName, et_vec(j), 'J2000', 'NONE', 'MOON');
    
    % Compute station-satellite vector in MCI
    rv_lander_sat_MCI = xx_MCI(j,:)' - rv_lander_MCI;
    
    % Transformation from ECI to topocentric frame
    ROT_MCI2TOPO = cspice_sxform('J2000', lander_TOPO, et_vec(j));

    % Convert state into topocentric frame
    rv_lander_sat_topo = ROT_MCI2TOPO*rv_lander_sat_MCI;

    % Compute range, azimuth and elevation using cspice_xfmsta
    rll_station_sat = cspice_xfmsta(rv_lander_sat_topo,'RECTANGULAR','LATITUDINAL','MOON');
    
    range(j)   = rll_station_sat(1);   % [km]
    azimuth(j) = rll_station_sat(2);   % [rad]
    elevation(j) = rll_station_sat(3); % [rad]
   
end
    if nargin == 4
       % Gaussian Error
       range = mvnrnd(range,sigma^2);
    end

end

%------------------------------------------------------------------------
function [xx, Px] = UKF(xx_0, P0, Weight, lambda, tspan, mu, measurments, sigma_mat, Mradii)
%----------------------------------------------------------------
% Unscented Kalman Filter (UKF), this function implements a sequential
% filter which propagate given initial state of N dimension with its relative covariance
% matrix and gives a progressivly better estimate of the state by
% introducing a provided set of m measurments
%
% INPUT:
% xx_0         [Nx1] Initial state vector
% P0           [NxN] Initial covariance matrix
% Weight       [struct] Containing the wheight matrices of the UT transformation
% lamda        [1] Lambda parameter of the UT
% tspan        [nx1] Time vector 
% mu           [1] Lunar gravity parameter
% measurments  [nxm] Measurment matrix
% sigma_mat    [mxm] Variance matrix of measurments noise
% Mradii       [1] Lunar equatorial radius [km]
%
% OUTPUT:
% xx           [nxN] State over time
% Px           [NxNxn] Covariance matrix over time
% -------------------------------------------------------------------------
    if nargin == 9
      Mradii_eq = Mradii;
    end

    W_mean = Weight.mean;
    W_covariance = Weight.covariance;
    P_old = P0;
    xx_old = xx_0';
    n = length(xx_0);
    xx = zeros(length(tspan),n);
    Px = zeros(n,n,length(tspan));
    y_propSigmaPoints = zeros(size(measurments,2), 2*n +1);

    xx(1,:) = xx_0;
    Px(:,:,1) = P0;
    
    for i = 1:length(tspan)-1

        % Pre-Processing for UT
        sigmaPoints = zeros(n,1+2*n);
        sqrt_Mat = chol((n+lambda)*P_old,"lower");
        propSigmaPoints =  zeros(n,1+2*n);
    
        % Sigma Points generation
        sigmaPoints(:,1) =  xx_old;
        sigmaPoints(:,2:n+1) = repmat(xx_old, 1, n) + sqrt_Mat;  
        sigmaPoints(:,n+2:2*n+1) = repmat(xx_old, 1, n) - sqrt_Mat; 

        % Sigma Points propagation
        options = odeset('RelTol',1e-13,'AbsTol',1e-13);
        time = [tspan(i) tspan(i+1)];

        for j = 1:1+2*n
            [~,sigma_prop] = ode113(@(t,y) ode_2bp(t,y,mu), time, sigmaPoints(1:6,j),options);
            propSigmaPoints(1:6,j) = sigma_prop(end,:)';

            if n > 6  % for ex 3.4
                % Longitude, Latitude
                propSigmaPoints(7:end,j) = sigmaPoints(7:end,j);
                
                % Transform sigma points related to Longitude and Latitude
                % into the measurment space (range)

                % Conversion (MCMF)
                rLand_IAU = cspice_latrec(Mradii_eq, propSigmaPoints(end-1,j), propSigmaPoints(end,j));
                % Rotation from MCMF to J2000 
                rotMat = cspice_sxform('IAU_MOON', 'J2000', tspan(i));
                % Lander coordinates in MCI from lon,lat
                xxLand_MCI = rotMat* [rLand_IAU; 0; 0; 0];
                % Transformation into measurment space
                [~, ~, range] = visibility(xxLand_MCI', propSigmaPoints(1:6,j)', tspan(i));
                
                y_propSigmaPoints(4,j) = range;
            end
        end
        
        % A priori mean and covariance
        xx_min = sum(W_mean.*propSigmaPoints,2);
        P_min = (propSigmaPoints - repmat(xx_min, 1, 2*n+1))*diag(W_covariance)*(propSigmaPoints - repmat(xx_min, 1, 2*n+1))';
       
        % Get sigma point related to measurments
        y_propSigmaPoints(1:3,:) = propSigmaPoints(1:3,:);
    
        % Mean of the measurments
        yMean_propSigmaPoints = sum(W_mean.*y_propSigmaPoints,2);
        
        % Update step
        P_yy = (y_propSigmaPoints - repmat(yMean_propSigmaPoints, 1, 2*n+1))*diag(W_covariance)*(y_propSigmaPoints - repmat(yMean_propSigmaPoints, 1, 2*n+1))'...
               + sigma_mat;
        P_xy = (propSigmaPoints - repmat(xx_min, 1, 2*n+1))*diag(W_covariance)*(y_propSigmaPoints - repmat(yMean_propSigmaPoints, 1, 2*n+1))';

        % Kalman Gain
        KK = (P_yy \ P_xy')';

        % A posteriori mean and covariance
        xx_plus = xx_min + KK*(measurments(i+1,:)' - yMean_propSigmaPoints);
        P_plus = P_min - KK*P_yy*KK';
  
        % Update State and Covariance
        xx(i+1,:) = xx_plus;
        Px(:,:,i+1) = P_plus;
    
        % Update variable
        P_old = P_plus;
        xx_old = xx_plus;
    end
    
    
end

%------------------------------------------------------------------------
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