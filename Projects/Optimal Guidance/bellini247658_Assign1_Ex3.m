% Spacecraft Guidance and Navigation (2024/2025)
% Assignment # 1, Exercise 3
% Author: Davide Bellini

clc; clearvars; close all; cspice_kclear()

% Load kernels
cspice_furnsh('ex02.tm');
% Plot settings
plot_set()

%% DATA
% NOTATION: Dimensional data with the index "_d"

h_i_d = 800;           % [km] Altitude of departure orbit
h_f_d = 1000;          % [km] Altitude of arrival orbit
del_i_d = 0.75;        % [deg] Inclination change
R_E_d = 6378.1366;     % [km] Earth radius
mu_d = 398600.435;     % [km^3/s^2] Earth gravitational parameter
rho_0_d = 750 + R_E_d; % [km] Reference radius for debris flux
k1 = 1e-5;             % [DU^−1] Debris spatial density constant 1
k2 = 1e-4;             % [DU^2] Debris spatial density constant 2
m_0_d = 1000;          % [kg] Initial mass
Tmax_d = 3;            % [N] Maximum thrust
Isp_d = 3120;          % [s] Specific impulse
DU = 7178.1366;        % [km] Distance Unit
MU = m_0_d;            % [kg] Mass Unit

% Debries spatial density
q_d = @(rho) k1./ (k2 + ((rho - rho_0_d)./DU).^2);
% TU, VU
TU = sqrt(DU^3/mu_d);
VU = DU/TU;
% Initial and final radius
r_i_d = R_E_d + h_i_d; %[km]
r_f_d = R_E_d + h_f_d; %[km]

%% --------------------------- 3.1 -----------------------
% Debries spatial density

% Plot debries spatial density
h_vec = linspace(r_i_d- 100, r_f_d + 100,10000);
Density = q_d(h_vec)*MU/DU^3;
figure
plot(h_vec- R_E_d ,Density,'LineWidth',2)
grid on
xlabel('Altitude [$km$]')
ylabel('Debries spatial density [$kg/km^3$]')
legend('Debris Density')

% Intial and final states J2000
% INITIAL
x_i = r_i_d; 
y_i = 0; 
z_i = 0;
vx_i = 0;
vy_i = sqrt(mu_d/r_i_d); 
vz_i = 0;
xx_i_dim = [x_i y_i z_i vx_i vy_i vz_i]'; 
% FINAL
x_f = r_f_d; 
y_f = 0; 
z_f = 0;
vx_f = 0;
vy_f = sqrt(mu_d/r_f_d)*cos( deg2rad(del_i_d));
vz_f = sqrt(mu_d/r_f_d)*sin( deg2rad(del_i_d));
xx_f_dim = [x_f y_f z_f vx_f vy_f vz_f]'; 

%% ----------------------------- 3.2 -----------------------------------
% Adimensionalization
adim_state = [DU DU DU VU VU VU]';
xx_i = xx_i_dim./adim_state; % initial state
xx_f = xx_f_dim./adim_state; % Final state

h_i = h_i_d/DU;                % [-] scaled altitude of departure orbit 
h_f = h_f_d/DU;                % [-] scaled altitude of arrival orbit
del_i = deg2rad(del_i_d);      % [-] scaled nclination change
R_E = R_E_d/DU;                % [-] scaled Earth radius
mu = 1;                        % [-] scaled Earth gravitational parameter
rho_0 = rho_0_d/DU;            % [-] scaled reference radius for debris flux
m_0 = m_0_d/MU;                % [-] Initial mass
Tmax = Tmax_d/(MU*DU*1000)*TU^2; % [-] Maximum thrust
Isp = Isp_d/TU;                               % [-] Specific impulse
g_0 = 9.81 *1e-3*TU^2/DU;                  % [-] gravity constant
q = @(rho) k1 ./ (k2 + (rho - rho_0).^2); % [-] debries spatial density
u = 1;

%% ----------------------------- 3.3 -------------------------------------
% Initial conditions
xx_0 = [xx_i;m_0];
% Target parameters
lam_m_f = 0;
H_f = 0;
target = [xx_f; lam_m_f; H_f];

% Intial condition while
n_iter = 1;
n_max = 15; % find the target solution around 20*pi could be required more the one attemp
condition = 0;
% time
tf_0 = 20*pi;
%  3*sigma distribution around target solution
sigma = (19.5*pi - 20.5*pi) / 6;
% Option
options = optimoptions('fsolve', 'MaxFunctionEvaluations', 10000, 'MaxIterations', 1000, 'Display','none',...
    'Algorithm', 'levenberg-marquardt', 'FunctionTolerance', 1e-10, 'OptimalityTolerance', 1e-08, 'StepTolerance', 1e-08,...
    'SpecifyObjectiveGradient',true,'UseParallel',false);

% FIND LAMBDA, TF
while n_iter <= n_max && condition == 0

    lambda_span = [-250 250];
    lambda_0 = zeros(7,1);
    lambda_0(1:6) = (lambda_span(2)-lambda_span(1)).*rand(6,1) +lambda_span(1);
    lambda_0(7) = lambda_span(2).*rand(1,1);
    tf_0 = tf_0 +rand(1,1)*sigma;
    ll_t_0 = [lambda_0;tf_0];

    % % To avoid long computaion the reader could use this previously found
    % % values by perturbing a valid solution.
    % ll_t_0 = [  -214.981221920134;
    % -10.3658727476080;
    % 0.885567778817514;
    % -10.3929207472605;
    % -214.610452407053;
    % -112.945358139111;
    % 2.59644920755792;
    % 64.4801062409445];
   

% zero finding problem
ObjFun = @(ll_t) Lamda_tf(ll_t,xx_0,mu,u,Isp,g_0,k1,k2,rho_0,Tmax,target);
%[exfg,err_grad] = checkGradients(ObjFun,ll_t_0,'Display','on');
%
[ll_t,lam_tf_err,exitflag,output] = fsolve(ObjFun,ll_t_0,options);
%count
n_iter = n_iter+1;

% Check the solution
if ll_t(end) > 21*pi || ll_t(end) < 19*pi
    condition = 0;
else
    condition = 1;
end

end

% HAMILTONIAN and PRIMER VECTOR
% Initial conditions
Phi0 = eye(14);
% Append to initial conditions the conditions for the STM
xx_ll_0 = [xx_0; ll_t(1:7); Phi0(:)];
% Propagation with STM
ode_options =  odeset('reltol', 1e-12, 'abstol', 1e-12);
[tt,results] = ode113(@(t,x) EL_eq_STM(t,x,mu,u,Isp,g_0,k1,k2,rho_0,Tmax), [0 ll_t(end)], xx_ll_0, ode_options);

% Post-processing
transf_duration = ll_t(end)*TU/60;              % transfer duration [min]
finalMass = results(end,7)*MU;                  % final mass [kg]
errPos = norm(xx_f(1:3) - results(end,1:3)')*DU; % Position error norm [km]
errVel = norm(xx_f(4:6) - results(end,4:6)')*VU*1e3; % Velocity error norm [m/s]
%
figure
plot3(results(:,1),results(:,2),results(:,3),'LineWidth',1.2)
hold on
plot3(results(1,1),results(1,2),results(1,3),'Marker','square','LineWidth',3,'LineStyle','none')
plot3(results(end,1),results(end,2),results(end,3),'Marker','square','LineWidth',3,'LineStyle','none')
grid on
xlabel('x [DU]')
ylabel('y [DU]')
zlabel('z [DU]')
title('Transfer orbit T = 3','FontSize',15)
legend('Transfer orbit','Starting point', 'Arrival point')

figure
plot3(results(:,1),results(:,2),results(:,3),'LineWidth',0.5)
hold on
plot3(results(1,1),results(1,2),results(1,3),'Marker','square','LineWidth',3,'LineStyle','none')
plot3(results(end,1),results(end,2),results(end,3),'Marker','square','LineWidth',3,'LineStyle','none')
grid on
xlabel('x [DU]')
ylabel('y [DU]')
zlabel('z [DU]')
title('Transfer orbit T = 3','FontSize',15)
legend('Transfer orbit','Starting point', 'Arrival point')
axis equal
% Hamiltonian 
H_time = NaN(length(tt),1);
for j = 1:length(tt)
    H_j = Hamiltonian(results(j,1:14),mu,u,Isp,g_0,k1,k2,rho_0,Tmax);
    H_time(j) = H_j;
end

figure()
plot(tt,H_time,'Color','b','LineWidth',1.5)
grid on
title('Hamiltonian over time')
ylabel('Hamiltonian [-]')
xlabel('Time [TU]')
legend('Hamiltonian')

% PRIME VECTOR   
alpha_NTW = NaN(length(tt),3);
for j = 1:length(tt)
    % Rotation
    lv_vec = results(j,11:13);
    alpha_ECI = -lv_vec./norm(lv_vec);
    alpha_NTW(j,:) = ECI2NTW(results(j,:),alpha_ECI);
end

figure()
plot(tt,alpha_NTW(:,1),'LineWidth',1.2)
hold on
plot(tt,alpha_NTW(:,2),'LineWidth',1.2)
plot(tt,alpha_NTW(:,3),'LineWidth',1.2)
legend('N','T','W')
title('Time evolution of the primer vector')
xlabel('Time [TU]')
ylabel('Primer direction [-]')
grid on

%% -------------------------- 3.4 ----------------------------------------
% Numerical continuation

% Adimensionalization
T_max_new = 2.86/(MU*DU*1000)*TU^2;
step = - 0.01/(MU*DU*1000)*TU^2;
T_span = Tmax : step : T_max_new;
ll_t_iter = NaN(8,length(T_span));
ll_t_new = ll_t; 

options = optimoptions('fsolve', 'MaxFunctionEvaluations', 10000, 'MaxIterations', 1000, 'Display','iter',...
    'Algorithm', 'trust-region-dogleg', 'FunctionTolerance', 1e-10, 'OptimalityTolerance', 1e-08, 'StepTolerance', 1e-08,...
    'SpecifyObjectiveGradient',true,'UseParallel',false);

for jj = 1:length(T_span)
    % Zero finding problem
    ObjFun = @(ll_t) Lamda_tf(ll_t,xx_0,mu,u,Isp,g_0,k1,k2,rho_0,T_span(jj),target);
    [ll_t_low] = fsolve(ObjFun,ll_t_new,options);
    % Save data
    ll_t_iter(:,n_iter) = ll_t_low;
    %
    ll_t_new = ll_t_low;

end

% HAMILTONIAN and PRIMER VECTOR
% Initial conditions
Phi0 = eye(14);
% Append to initial conditions the conditions for the STM
xx_ll_0_new = [xx_0; ll_t_new(1:7); Phi0(:)];
% Propagation with STM
ode_options =  odeset('reltol', 1e-12, 'abstol', 1e-12);
[tt_new,results_new] = ode113(@(t,x) EL_eq_STM(t,x,mu,u,Isp,g_0,k1,k2,rho_0,T_max_new), [0 ll_t_new(end)], xx_ll_0_new, ode_options);

% Post-processing
transf_duration_new = ll_t_new(end)*TU/60;
finalMass_new = results_new(end,7)*MU;
errPos_new = norm(xx_f(1:3) - results_new(end,1:3)')*DU;
errVel_new = norm(xx_f(4:6) - results_new(end,4:6)')*VU*1e3;

%
figure
plot3(results_new(:,1),results_new(:,2),results_new(:,3),'LineWidth',1.2)
hold on
plot3(results_new(1,1),results_new(1,2),results_new(1,3),'Marker','square','LineWidth',3,'LineStyle','none')
plot3(results_new(end,1),results_new(end,2),results_new(end,3),'Marker','square','LineWidth',3,'LineStyle','none')
grid on
xlabel('x [DU]')
ylabel('y [DU]')
zlabel('z [DU]')
title('Transfer orbit T = 2.86','FontSize',15)
legend('Transfer orbit T = 2.86','Starting point', 'Arrival point')

% Hamiltonian 
H_time_new = NaN(length(tt_new),1);
for j = 1:length(tt_new)
    H_j = Hamiltonian(results_new(j,1:14),mu,u,Isp,g_0,k1,k2,rho_0,T_max_new);
    H_time_new(j) = H_j;
end

figure()
plot(tt_new,H_time_new,'Color','b','LineWidth',1.5)
grid on
title('Hamiltonian over time after numerical continuation')
ylabel('Hamiltonian [-]')
xlabel('Time [TU]')
legend('Hamiltonian')

% PRIME VECTOR   
alpha_NTW_new = NaN(length(tt_new),3);
for j = 1:length(tt_new)
    % Rotation
    lv_vec = results_new(j,11:13);
    alpha_ECI = -lv_vec./norm(lv_vec);
    alpha_NTW_new(j,:) = ECI2NTW(results_new(j,:),alpha_ECI);
end

figure()
plot(tt_new,alpha_NTW_new(:,1),'LineWidth',1.2)
hold on
plot(tt_new,alpha_NTW_new(:,2),'LineWidth',1.2)
plot(tt_new,alpha_NTW_new(:,3),'LineWidth',1.2)
legend('N','T','W')
title('Time evolution of the primer vector T = 2.86')
xlabel('Time [TU]')
ylabel('Primer direction [-]')
grid on


cspice_kclear()


%% FUNCTION
% -----------------------------------------------------------------
function H = Hamiltonian(xx_ll_t,mu,u,Isp,g_0,k1,k2,rho_0,Tmax)
%------------------------------------------------------------------------
% Computes Hamiltonian of the problem 
%
% INPUT
% xx_ll_t   [14x1] State, Costate at a given time
% mu        [1] Gravitational parameter of the attractor
% u         [1] Thrust throttle level
% Isp       [1] Specific impulse
% g_0       [1] Gravity costant at sea level
% k1        [1] Debries density parameter
% k2        [1] Debries density parameter
% rho_0     [1] Reference radius
% Tmax      [1] Thrust
%
% OUTPUT
% H         [1] Hamiltonian
% ------------------------------------------------------------------------
% Notation
rr = xx_ll_t(1:3)';
vv = xx_ll_t(4:6)';
m = xx_ll_t(7);
lr = xx_ll_t(8:10)';
lv = xx_ll_t(11:13)';
lm = xx_ll_t(14)';
r = norm(rr);
% primer vector
alpha = -lv/norm(lv);
lam_v = [lr ;lv ;lm];
% EOM
eom = [vv; -mu/r^3 *rr + u*Tmax/m*alpha; -u*Tmax/(Isp*g_0)];
% Debries spatial density
q =  k1 / (k2 + (r - rho_0)^2);
% Hamiltonian
H = q + lam_v'*eom;
end

% -----------------------------------------------------------------
function xx_ll = EL_eq_STM(~,xx_ll,mu,u,Isp,g_0,k1,k2,rho_0,Tmax)
%------------------------------------------------------------------------
% Computes the Euler Lagrange equation for the PMP, together with the STM
% of state and costate
%
% INPUT
% xx_ll     [210x1] State, Costate and STM
% mu        [1] Gravitational parameter of the attractor
% u         [1] Thrust throttle level
% Isp       [1] Specific impulse
% g_0       [1] Gravity costant at sea level
% k1        [1] Debries density parameter
% k2        [1] Debries density parameter
% rho_0     [1] Reference radius
% Tmax      [1] Thrust
%
% OUTPUT
% lam_tf_err   [210x1] RHS of the EL equation and STM
% ------------------------------------------------------------------------

% Notation
rr = xx_ll(1:3);
r = norm(rr);
vv = xx_ll(4:6);
m = xx_ll(7);
lr = xx_ll(8:10);
lv = xx_ll(11:13);

% primer vector
alpha = -lv/norm(lv);
% derivative of the debries spatial density
dq = [(2*k1*rr(1)*(rho_0 - (rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2)))/((k2 + (rho_0 - (rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2))^2)^2*(rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2));
      (2*k1*rr(2)*(rho_0 - (rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2)))/((k2 + (rho_0 - (rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2))^2)^2*(rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2));
      (2*k1*rr(3)*(rho_0 - (rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2)))/((k2 + (rho_0 - (rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2))^2)^2*(rr(1)^2 + rr(2)^2 + rr(3)^2)^(1/2))];

x = rr(1); y = rr(2); z = rr(3); lvx = lv(1); lvy = lv(2); lvz = lv(3);
% Put PHI in matrix form
PHI = reshape(xx_ll(15:end),14,14);
% dfdx
dfdx = [                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            0, 1, 0, 0,                                              0,  0,  0,  0,                                                         0,                                             0,                                                         0, 0;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            0, 0, 1, 0,                                              0,  0,  0,  0,                                                         0,                                                         0,                                                         0, 0;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            0, 0, 0, 1,                                              0,  0,  0,  0,                                                         0,                                                         0,                                                         0, 0;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         -(mu*(- 2*x^2 + y^2 + z^2))/(x^2 + y^2 + z^2)^(5/2),                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           (3*mu*x*y)/(x^2 + y^2 + z^2)^(5/2),                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           (3*mu*x*z)/(x^2 + y^2 + z^2)^(5/2), 0, 0, 0, (Tmax*lvx)/(m^2*(lvx^2 + lvy^2 + lvz^2)^(1/2)),  0,  0,  0, -(Tmax*(lvy^2 + lvz^2))/(m*(lvx^2 + lvy^2 + lvz^2)^(3/2)),          (Tmax*lvx*lvy)/(m*(lvx^2 + lvy^2 + lvz^2)^(3/2)),          (Tmax*lvx*lvz)/(m*(lvx^2 + lvy^2 + lvz^2)^(3/2)), 0;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          (3*mu*x*y)/(x^2 + y^2 + z^2)^(5/2),                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            -(mu*(x^2 - 2*y^2 + z^2))/(x^2 + y^2 + z^2)^(5/2),                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           (3*mu*y*z)/(x^2 + y^2 + z^2)^(5/2), 0, 0, 0, (Tmax*lvy)/(m^2*(lvx^2 + lvy^2 + lvz^2)^(1/2)),  0,  0,  0,          (Tmax*lvx*lvy)/(m*(lvx^2 + lvy^2 + lvz^2)^(3/2)), -(Tmax*(lvx^2 + lvz^2))/(m*(lvx^2 + lvy^2 + lvz^2)^(3/2)),          (Tmax*lvy*lvz)/(m*(lvx^2 + lvy^2 + lvz^2)^(3/2)), 0;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          (3*mu*x*z)/(x^2 + y^2 + z^2)^(5/2),                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           (3*mu*y*z)/(x^2 + y^2 + z^2)^(5/2),                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            -(mu*(x^2 + y^2 - 2*z^2))/(x^2 + y^2 + z^2)^(5/2), 0, 0, 0, (Tmax*lvz)/(m^2*(lvx^2 + lvy^2 + lvz^2)^(1/2)),  0,  0,  0,          (Tmax*lvx*lvz)/(m*(lvx^2 + lvy^2 + lvz^2)^(3/2)),          (Tmax*lvy*lvz)/(m*(lvx^2 + lvy^2 + lvz^2)^(3/2)), -(Tmax*(lvx^2 + lvy^2))/(m*(lvx^2 + lvy^2 + lvz^2)^(3/2)), 0;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            0, 0, 0, 0,                                              0,  0,  0,  0,                                                         0,                                                         0,                                                         0, 0;
(15*mu*x^2*(lvx*x + lvy*y + lvz*z))/(x^2 + y^2 + z^2)^(7/2) - (3*mu*(lvx*x + lvy*y + lvz*z))/(x^2 + y^2 + z^2)^(5/2) - (2*k1*(rho_0 - (x^2 + y^2 + z^2)^(1/2)))/((k2 + (rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)^(1/2)) + (2*k1*x^2)/((k2 + (rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)) - (6*lvx*mu*x)/(x^2 + y^2 + z^2)^(5/2) + (2*k1*x^2*(rho_0 - (x^2 + y^2 + z^2)^(1/2)))/((k2 + (rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)^(3/2)) - (8*k1*x^2*(rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)/((k2 + (rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)^3*(x^2 + y^2 + z^2)),                                                                                                                                   (2*k1*x*y)/((k2 + (rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)) - (3*lvx*mu*y)/(x^2 + y^2 + z^2)^(5/2) - (3*lvy*mu*x)/(x^2 + y^2 + z^2)^(5/2) + (15*mu*x*y*(lvx*x + lvy*y + lvz*z))/(x^2 + y^2 + z^2)^(7/2) + (2*k1*x*y*(rho_0 - (x^2 + y^2 + z^2)^(1/2)))/((k2 + (rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)^(3/2)) - (8*k1*x*y*(rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)/((k2 + (rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)^3*(x^2 + y^2 + z^2)),                                                                                                                                   (2*k1*x*z)/((k2 + (rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)) - (3*lvx*mu*z)/(x^2 + y^2 + z^2)^(5/2) - (3*lvz*mu*x)/(x^2 + y^2 + z^2)^(5/2) + (15*mu*x*z*(lvx*x + lvy*y + lvz*z))/(x^2 + y^2 + z^2)^(7/2) + (2*k1*x*z*(rho_0 - (x^2 + y^2 + z^2)^(1/2)))/((k2 + (rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)^(3/2)) - (8*k1*x*z*(rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)/((k2 + (rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)^3*(x^2 + y^2 + z^2)), 0, 0, 0,                                              0,  0,  0,  0,        (mu*(- 2*x^2 + y^2 + z^2))/(x^2 + y^2 + z^2)^(5/2),                       -(3*mu*x*y)/(x^2 + y^2 + z^2)^(5/2),                       -(3*mu*x*z)/(x^2 + y^2 + z^2)^(5/2), 0;
                                                                                                                                  (2*k1*x*y)/((k2 + (rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)) - (3*lvx*mu*y)/(x^2 + y^2 + z^2)^(5/2) - (3*lvy*mu*x)/(x^2 + y^2 + z^2)^(5/2) + (15*mu*x*y*(lvx*x + lvy*y + lvz*z))/(x^2 + y^2 + z^2)^(7/2) + (2*k1*x*y*(rho_0 - (x^2 + y^2 + z^2)^(1/2)))/((k2 + (rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)^(3/2)) - (8*k1*x*y*(rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)/((k2 + (rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)^3*(x^2 + y^2 + z^2)), (15*mu*y^2*(lvx*x + lvy*y + lvz*z))/(x^2 + y^2 + z^2)^(7/2) - (3*mu*(lvx*x + lvy*y + lvz*z))/(x^2 + y^2 + z^2)^(5/2) - (2*k1*(rho_0 - (x^2 + y^2 + z^2)^(1/2)))/((k2 + (rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)^(1/2)) + (2*k1*y^2)/((k2 + (rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)) - (6*lvy*mu*y)/(x^2 + y^2 + z^2)^(5/2) + (2*k1*y^2*(rho_0 - (x^2 + y^2 + z^2)^(1/2)))/((k2 + (rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)^(3/2)) - (8*k1*y^2*(rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)/((k2 + (rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)^3*(x^2 + y^2 + z^2)),                                                                                                                                   (2*k1*y*z)/((k2 + (rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)) - (3*lvy*mu*z)/(x^2 + y^2 + z^2)^(5/2) - (3*lvz*mu*y)/(x^2 + y^2 + z^2)^(5/2) + (15*mu*y*z*(lvx*x + lvy*y + lvz*z))/(x^2 + y^2 + z^2)^(7/2) + (2*k1*y*z*(rho_0 - (x^2 + y^2 + z^2)^(1/2)))/((k2 + (rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)^(3/2)) - (8*k1*y*z*(rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)/((k2 + (rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)^3*(x^2 + y^2 + z^2)), 0, 0, 0,                                              0,  0,  0,  0,                       -(3*mu*x*y)/(x^2 + y^2 + z^2)^(5/2),          (mu*(x^2 - 2*y^2 + z^2))/(x^2 + y^2 + z^2)^(5/2),                       -(3*mu*y*z)/(x^2 + y^2 + z^2)^(5/2), 0;
                                                                                                                                  (2*k1*x*z)/((k2 + (rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)) - (3*lvx*mu*z)/(x^2 + y^2 + z^2)^(5/2) - (3*lvz*mu*x)/(x^2 + y^2 + z^2)^(5/2) + (15*mu*x*z*(lvx*x + lvy*y + lvz*z))/(x^2 + y^2 + z^2)^(7/2) + (2*k1*x*z*(rho_0 - (x^2 + y^2 + z^2)^(1/2)))/((k2 + (rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)^(3/2)) - (8*k1*x*z*(rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)/((k2 + (rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)^3*(x^2 + y^2 + z^2)),                                                                                                                                   (2*k1*y*z)/((k2 + (rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)) - (3*lvy*mu*z)/(x^2 + y^2 + z^2)^(5/2) - (3*lvz*mu*y)/(x^2 + y^2 + z^2)^(5/2) + (15*mu*y*z*(lvx*x + lvy*y + lvz*z))/(x^2 + y^2 + z^2)^(7/2) + (2*k1*y*z*(rho_0 - (x^2 + y^2 + z^2)^(1/2)))/((k2 + (rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)^(3/2)) - (8*k1*y*z*(rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)/((k2 + (rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)^3*(x^2 + y^2 + z^2)), (15*mu*z^2*(lvx*x + lvy*y + lvz*z))/(x^2 + y^2 + z^2)^(7/2) - (3*mu*(lvx*x + lvy*y + lvz*z))/(x^2 + y^2 + z^2)^(5/2) - (2*k1*(rho_0 - (x^2 + y^2 + z^2)^(1/2)))/((k2 + (rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)^(1/2)) + (2*k1*z^2)/((k2 + (rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)) - (6*lvz*mu*z)/(x^2 + y^2 + z^2)^(5/2) + (2*k1*z^2*(rho_0 - (x^2 + y^2 + z^2)^(1/2)))/((k2 + (rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)^2*(x^2 + y^2 + z^2)^(3/2)) - (8*k1*z^2*(rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)/((k2 + (rho_0 - (x^2 + y^2 + z^2)^(1/2))^2)^3*(x^2 + y^2 + z^2)), 0, 0, 0,                                              0,  0,  0,  0,                       -(3*mu*x*z)/(x^2 + y^2 + z^2)^(5/2),                       -(3*mu*y*z)/(x^2 + y^2 + z^2)^(5/2),          (mu*(x^2 + y^2 - 2*z^2))/(x^2 + y^2 + z^2)^(5/2), 0;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            0, 0, 0, 0,                                              0, -1,  0,  0,                                                         0,                                                         0,                                                         0, 0;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            0, 0, 0, 0,                                              0,  0, -1,  0,                                                         0,                                                         0,                                                         0, 0;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            0, 0, 0, 0,                                              0,  0,  0, -1,                                                         0,                                                         0,                                                         0, 0;

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            0, 0, 0, 0,     (2*Tmax*(lvx^2 + lvy^2 + lvz^2)^(1/2))/m^3,  0,  0,  0,           -(Tmax*lvx)/(m^2*(lvx^2 + lvy^2 + lvz^2)^(1/2)),           -(Tmax*lvy)/(m^2*(lvx^2 + lvy^2 + lvz^2)^(1/2)),           -(Tmax*lvz)/(m^2*(lvx^2 + lvy^2 + lvz^2)^(1/2)), 0];

% Compute the derivative of the STM
PHIdot = dfdx*PHI;
xx_ll = zeros(210,1);
% EOM
xx_ll(1:3)   = vv;
xx_ll(4:6)   = - mu/r^3 * rr + u*Tmax/m * alpha;
xx_ll(7)     = -u*Tmax/(Isp*g_0);
xx_ll(8:10)  = -dq + mu/r^3*lv - 3*mu/r^5 *(dot(rr,lv))*rr;
xx_ll(11:13) = -lr;
xx_ll(14)    = -u*Tmax/m^2*norm(lv);
xx_ll(15:end) = PHIdot(:);
end

%---------------------------------------------------------------------
function [lam_tf_err,fsolve_grad] = Lamda_tf(ll_t,xx_0,mu,u,Isp,g_0,k1,k2,rho_0,Tmax,target)
%------------------------------------------------------------------------
% Objective function for the initial costate and final time optimization,
% Evaluate the residual from the rarget conditions.
%
% INPUT
% ll_t      [8x1] Costate and final time
% mu        [1] Gravitational parameter of the attractor
% u         [1] Thrust throttle level
% Isp       [1] Specific impulse
% g_0       [1] Gravity costant at sea level
% k1        [1] Debries density parameter
% k2        [1] Debries density parameter
% rho_0     [1] Reference radius
% Tmax      [1] Thrust
% Target    [8x1] Target conditions
%
% OUTPUT
% lam_tf_err   [8x1] Residuals
% fsolve_grad  [8x8] Derivatives of the objective function
% ------------------------------------------------------------------------

% Options
ode_options =  odeset('reltol', 1e-12, 'abstol', 1e-12);
% Initialization
t_end = ll_t(end);
xx_ll_in = [xx_0;ll_t(1:7)];
Phi0 = eye(14);

% Append to initial conditions the conditions for the STM
xx_ll_0 = [xx_ll_in; Phi0(:)];  
% propagation
[~,xx_ll_v] = ode113(@(t,x) EL_eq_STM(t,x,mu,u,Isp,g_0,k1,k2,rho_0,Tmax), [0 t_end], xx_ll_0, ode_options);
xx_ll_f = xx_ll_v(end,1:14);
PHI = reshape(xx_ll_v(end,15:end),14,14);
% Hamiltonian
H = Hamiltonian(xx_ll_f,mu,u,Isp,g_0,k1,k2,rho_0,Tmax);
% Gradient
fsolve_grad = Gradient(PHI,xx_ll_f,mu,u,Isp,g_0,k1,k2,rho_0,Tmax);
% errors
lam_tf_err =  [xx_ll_f(1:6)'; xx_ll_f(end); H]- target ;
end

%---------------------------------------------------------------------
function fsolve_grad = Gradient(PHI,xx_ll_f,mu,u,Isp,g_0,k1,k2,rho_0,Tmax)
%------------------------------------------------------------------------
% Compute the derivatives of the objective function with respect to the
% state variables
%
% INPUT
% PHI       [14x14] State transition matrix of state and costate
% xx_ll_f   [14x1] state and costate
% mu        [1] Gravitational parameter of the attractor
% u         [1] Thrust throttle level
% Isp       [1] Specific impulse
% g_0       [1] Gravity costant at sea level
% k1        [1] Debries density parameter
% k2        [1] Debries density parameter
% rho_0     [1] Reference radius
% Tmax      [1] Thrust
%
% OUTPUT
% fsolve_grad  [8x8] Derivatives of the objective function
% ------------------------------------------------------------------------
       % Derivatives
       Gradient = zeros(8,8);
       xx_ll = [xx_ll_f';PHI(:)];
       xx_ll = EL_eq_STM(0,xx_ll,mu,u,Isp,g_0,k1,k2,rho_0,Tmax);

       Gradient(1:3,1:7) = PHI(1:3,8:end);
       Gradient(4:6,1:7) = PHI(4:6,8:end);
       Gradient(7,1:7)   = PHI(14,8:end);
       Gradient(1:8,end) = [xx_ll(1:6);xx_ll(14);0];
       Gradient(8,1:end-1) =  -xx_ll(8:14)'*PHI(1:7,8:14) + xx_ll(1:7)'*PHI(8:14,8:14) ;

       fsolve_grad = Gradient;
end

%---------------------------------------------------------------------
function NTW_vec = ECI2NTW(results,alpha_ECI)
%------------------------------------------------------------------------
% Rotate from ECI to NTW
%
% INPUT
% results   [nx6] State vector in ECI  
% alpha_ECI [nx3] Primer vector in ECI 
%
% OUTPUT
% NTW_vec   [nx3] Primer vector in NTW 
% ------------------------------------------------------------------------
% Get Coordinates
xx_ECI = results(1:3);
vv_ECI = results(4:6);
% Define the vectors
T = vv_ECI./norm(vv_ECI);
N = cross(xx_ECI,vv_ECI)./norm(cross(xx_ECI,vv_ECI));
W = cross(N,T);
% Rotation matrix
R = [N;T;W];
% Perform rotation
NTW_vec = R*alpha_ECI';
end

%-----------------------------------------------------------------------
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
