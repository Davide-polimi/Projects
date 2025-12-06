% Spacecraft Guidance and Navigation (2024/2025)
% Assignment # 1, Exercise 2
% Author: Davide Bellini

clearvars; close all; clc;
cspice_kclear()

% Load kernels
cspice_furnsh('ex02.tm');
% Plot settings
plot_set()
%% DATA
% Data from Table 4 
m_sun = 3.28900541*1e5;  % [-] Sun scaled mass
rho = 3.88811143*1e2;    % [-] Scaled Sun-(Earth-Moon) Batycenter distance
om_s = -9.25195985*1e-1; % [-] Scaled ancugal velocity of the sun
om_em = 2.66186135*1e-6; % [s^-1] Earth-Moon angular velocity
l_em = 3.84405*1e8;      % [m] Earth-Moon distance
h_i = 167;               % [km] Altitude of departure orbit
h_f = 100;               % [km] Altitude of arrival orbit

DU = 3.84405000*1e5;     % [km] Distance Unit
TU = 4.34256461;         % [days] Time Unit
VU = 1.02454018;         % [km/s] Velocity Unit

% First guess solution data
alpha = 0.2*pi;          % [rad] Location of the initial point on the circular orbit
beta = 1.41;             % [-] Initial to circular velocity ratio
t_i = 2;                 % [-] Initial time
del = 4;                 % [days] Transfer duration

%% --------------------- 2.1 -------------------------------
% First guess generation

% Computing mu
GM_E = cspice_bodvrd('Earth','GM',1);   % GM Earth
GM_M = cspice_bodvrd('Moon','GM',1);   % GM Moon
mu =  GM_M / (GM_E + GM_M);

% Extract Radii
E_RADII = cspice_bodvrd('Earth','RADII',3);
M_RADII = cspice_bodvrd('Moon','RADII',3);

% Create scaled variables   % Select equatorial radii
% Intitial conditions
r_i = (E_RADII(1) + h_i ) / DU;   % Scaled initial radius
v_circ_i = sqrt( (1-mu)/r_i);     % Scaled initial orbit circular velocity
% Final conditions
r_f = ( M_RADII(1) + h_f ) / DU;   % Scaled final radius
v_circ_f = sqrt(mu/r_f);           % Scaled final orbit circular velocity
% Create scaled earth moon barycenter velocity
om_scaled = om_em*l_em/ (VU*1000);

% FIRST GUESS
% first guess solution from the given parameters alpha, beta, ti, delta)
% find v_0
v_0 = beta*sqrt( (1-mu)/r_i);    
% find initial state x_0 , initial time t_i , final time t_f
t_i;
t_f = t_i + del;
x_0 = r_i*cos(alpha) - mu;
y_0 = r_i*sin(alpha);
vx_0 = -sin(alpha)*(v_0-r_i);
vy_0 = (v_0-r_i)*cos(alpha);
xx0 = [x_0, y_0, vx_0, vy_0]';

% Propagation
[~,~, xx, tt]  = propagate(t_i,xx0,t_f,mu,m_sun,rho,om_s);

% Earth position
P1 = [-mu; 0];
% Moon position
P2 = [1-mu, 0];

% Orbit recnosruction
theta = linspace(0,2*pi,3600); % One orbit
% Arrival and intial orbits
xx_orb_E = [-mu + r_i.*cos(theta); r_i.*sin(theta)];
xx_orb_M = [1 - mu + r_f.* cos(theta); r_f.*sin(theta)];

% FIGURE EM-barycenter RF
figure
plot(xx(:,1),xx(:,2),'LineWidth',1.5) % Plot the path of the initial guess
hold on
plot(xx_orb_E(1,:),xx_orb_E(2,:),'LineStyle','-.','Color',[0.4 0.4 0.4])
plot(xx_orb_M(1,:),xx_orb_M(2,:),'LineStyle','-.','Color',[0.4 0.4 0.4])
grid on
axis equal
plot(P1(1),P1(2),'Marker','o','LineWidth',1,'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','b')
plot(P2(1),P2(2),'Marker','o','LineWidth',1,'Color',[0.7, 0.7, 0.7],'MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor',[0.7, 0.7, 0.7])
xlabel('x [DU]')
ylabel('y [DU]')
title('Initial guess propagation in EM-barycenter','FontSize',15)

% Rotation to ECI
XX = EM2ECI(tt, xx, mu);
xx_ECI_E = [r_i.*cos(theta); r_i.*sin(theta)];
xx_Moon = [cos(theta); sin(theta)];

% FIGURE ECI
figure
plot(XX(:,1),XX(:,2),'LineWidth',1.5) % Plot the path of the initial guess
hold on
% Plot Earth
plot(0,0,'Marker','o','LineWidth',1,'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','b')
plot(xx_Moon(1,:),xx_Moon(2,:),'LineStyle','--','Color',[0.4 0.4 0.4])
grid on
xlabel('x [DU]')
ylabel('y [DU]')
title('Initial guess propagation in ECI','FontSize',15)
axis equal

%% ------------------------------ 2.2 -------------------------- 
% SIMPLE SHOOTING

% first guess
xx0_t = [xx0;t_i;t_f];

%--------------------------- Case a --------------------------------

% fmincon options
options = optimoptions('fmincon','Display','iter', 'Algorithm','active-set','ConstraintTolerance',1e-10,...
                   'OptimalityTolerance',1e-10); 
% Objective function
Obj_fun = @(x) Delta_V(x,mu,m_sun,rho,om_s, v_circ_i, v_circ_f);
% Constraints function
c = @(x) constraints(x,mu,m_sun,rho,om_s,r_i,r_f);

%  Without providing any derivative
[opt_state,DV,~,~,~,~,~] = fmincon(Obj_fun, xx0_t,[0 0 0 0 1 -1],0,[],[],[],[], c, options);
% Propagation
[~,~, xx, tt]  = propagate(opt_state(5),opt_state(1:4),opt_state(6),mu,m_sun,rho,om_s);

% FIGURE EM barycenter RF
figure
plot(xx(:,1),xx(:,2),'LineWidth',1.5) % Plot the path of the initial guess
hold on
plot(xx_orb_E(1,:),xx_orb_E(2,:),'LineStyle','-.','Color',[0.4 0.4 0.4])
plot(xx_orb_M(1,:),xx_orb_M(2,:),'LineStyle','-.','Color',[0.4 0.4 0.4])
grid on
axis equal
% Plot Earth and Moon
plot(P1(1),P1(2),'Marker','o','LineWidth',1,'MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','b')
plot(P2(1),P2(2),'Marker','o','LineWidth',1,'Color','b','MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','b')
xlabel('x [DU]')
ylabel('y [DU]')
title('Simple shooting propagation in EM','FontSize',15)
axis equal

% Figure ECI RF 
XX = EM2ECI(tt, xx, mu);
figure
plot(XX(:,1),XX(:,2),'LineWidth',1.5) % Plot the path of the initial guess
hold on
plot(xx_ECI_E(1,:),xx_ECI_E(2,:),'LineStyle','-.','Color',[0.4 0.4 0.4])
% Plot Earth
plot(0,0,'Marker','o','LineWidth',1,'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','b')
plot(xx_Moon(1,:),xx_Moon(2,:),'LineStyle','--','Color',[0.4 0.4 0.4])
grid on
xlabel('x [DU]')
ylabel('y [DU]')
title('Simple shooting propagation in ECI','FontSize',15)
axis equal


% ----------------------- Case b -------------------------------------
% Providing the derivatives to the solver

% Initial guess
xx0_t = [xx0;t_i;t_f];

% fmincon options
options = optimoptions('fmincon','Display','iter', 'Algorithm','active-set','ConstraintTolerance',1e-10,...
                       'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);
% Objective function
Obj_fun = @(x) Delta_V_STM(x,mu,m_sun,rho,om_s, v_circ_i, v_circ_f);
% [val, error] = checkGradients(Obj_fun,xx0_t,'Display','on'); % Used to control the provided gradient

% Constraints function
c = @(x) constraints_STM(x,mu,m_sun,rho,om_s,r_i,r_f);
%[val, error1] = checkGradients(c,xx0_t,IsConstraint=true);  % Used to control the provided gradient

% fmincon
[opt_state_s,DV_s] = fmincon(Obj_fun, xx0_t,[0 0 0 0 1 -1],0,[],[],[],[], c, options);

% Propagation
[xf,tf, xx, tt]  = propagate(opt_state_s(5),opt_state_s(1:4),opt_state_s(6),mu,m_sun,rho,om_s);

% FIGURE EM barycenter RF
figure
plot(xx(:,1),xx(:,2),'LineWidth',1.5) % Plot the path of the initial guess
hold on
plot(xx_orb_E(1,:),xx_orb_E(2,:),'LineStyle','-.','Color',[0.4 0.4 0.4])
plot(xx_orb_M(1,:),xx_orb_M(2,:),'LineStyle','-.','Color',[0.4 0.4 0.4])
grid on
axis equal
% Plot Earth and Moon
plot(P1(1),P1(2),'Marker','o','LineWidth',1,'MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','b')
plot(P2(1),P2(2),'Marker','o','LineWidth',1,'Color','b','MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','b')
xlabel('x [DU]')
ylabel('y [DU]')
title('Simple shooting with analytical derivatives propagation in EM, ','FontSize',15)

% Figure ECI RF 
XX = EM2ECI(tt, xx, mu);
figure
plot(XX(:,1),XX(:,2),'LineWidth',1.5) 
hold on
plot(xx_ECI_E(1,:),xx_ECI_E(2,:),'LineStyle','-.','Color',[0.4 0.4 0.4])
% Plot Earth
plot(0,0,'Marker','o','LineWidth',1,'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','b')
plot(xx_Moon(1,:),xx_Moon(2,:),'LineStyle','--','Color',[0.4 0.4 0.4])
grid on
xlabel('x [DU]')
ylabel('y [DU]')
title('Simple shooting with analytical derivatives in ECI','FontSize',15)
axis equal


%% ------------------------------- 2.3 ---------------------------------
% MULTIPLE SHOOTING
% first guess
N = 4;                    % Number of points
h = (t_f - t_i) / (N-1);  
t_vec = t_i:h:t_f;        % time vector
xx0_N = NaN(4,N-1);     
xx_jf = xx0;
% Find initial solutions
for j = 1 : length(t_vec)-1
    [xx_jf,tf_j, xx_j, tt_j]  = propagate(t_vec(j),xx_jf,t_vec(j+1),mu,m_sun,rho,om_s);
    xx0_N(:,j) = xx_jf;
end
% Initial guess multiple shooting
xx0_t = [xx0; xx0_N(:);t_i;t_f];
%
% fmincon options 
options = optimoptions('fmincon','Display','iter', 'Algorithm','active-set','ConstraintTolerance',1e-7,...
                       'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
                       'MaxIterations',1000,'MaxFunctionEvaluations',5000);
% Objective function
Obj_fun = @(x) Delta_V_ms(x,mu,v_circ_i, v_circ_f,N);
% [val, error] = checkGradients(Obj_fun,xx0_t,'Display','on'); % Used to check the provided gradient

% Constraints function
c = @(x) constraints_ms(x,mu,m_sun,rho,om_s,r_i,r_f,N,E_RADII/DU,M_RADII/DU);
% [val1, error1] = checkGradients(c,xx0_t,IsConstraint=true); % Used to check the provided gradient

% Multiple shooting
[opt_state_ms,DV_ms,exitflag,output,lambda,grad,hess_ms] = fmincon(Obj_fun, xx0_t,[zeros(1,4*N) 1 -1],0,[],[],[],[], c, options);

% Propagation
h = (opt_state_ms(end) - opt_state_ms(end-1)) / (N-1);
t_vec = opt_state_ms(end-1):h:opt_state_ms(end);
col = parula(N);

% Plot
figure
hold on
for jj = 1:N
   ind = (1:4)+ 4*(jj-1);
   if jj < N
   [xp_jj,tp_jj, xxp_jj, ttp_jj]  = propagate(t_vec(jj),opt_state_ms(ind),t_vec(jj+1),mu,m_sun,rho,om_s);
   plot(xxp_jj(:,1),xxp_jj(:,2),'Color',col(jj,:),'LineWidth',1.5) % Plot the path of the initial guess
   end
  k = 1 + 4*(jj-1);
  plot(opt_state_ms(k),opt_state_ms(k+1),'o','Marker','o','Color','r','LineWidth',2)
end
% Plot arrival and initial orbits
plot(xx_orb_E(1,:),xx_orb_E(2,:),'LineStyle','-.','Color',[0.4 0.4 0.4])
plot(xx_orb_M(1,:),xx_orb_M(2,:),'LineStyle','-.','Color',[0.4 0.4 0.4])
grid on
axis equal
% Plot Earth and Moon
plot(P1(1),P1(2),'Marker','o','LineWidth',1,'MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','b')
plot(P2(1),P2(2),'Marker','o','LineWidth',1,'Color','b','MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','b')
xlabel('x [DU]')
ylabel('y [DU]')
title('Multiple shooting with analytical derivatives in EM','FontSize',15)

%% ----------------------------------- 2.4 ---------------------------- 
% Theta finder and N-body propagator

% FIND THETA
theta_i = wrapTo2Pi(om_s*opt_state_s(end-1));
% EARTH-MOON Period
T_EM = 2*pi/om_em;
% STARTING EPOCH
epoch_0_str = 'September 28 00:00:00.000 TDB 2024';
et_0 = cspice_str2et(epoch_0_str);

% MAX EPOCH
et_max = et_0 + T_EM;
% Find the bound for the zero finding problem 
t_span = linspace(et_0,et_max,T_EM);
theta_fun = @(t) Thetafinder(t,theta_i);
theta_vec = theta_fun(t_span);

% Plot to find the bounds for the zero finding problem
% this plot is used only to adjust the search interval
figure
plot(theta_vec,'LineWidth',1.5)
grid on
xlabel('Time [s]')
ylabel('$\theta$ [rad]')
title('Evolution of the Sun-Earth-Moon angle')

% Zero finding problem
dt = 0.5*1e6; % arbitrary value from the previous plot to reduce the interval
t_bounds = [et_0 + dt, et_max];
options = optimset('Display', 'iter', 'TolFun', 1e-12, 'TolX', 1e-12);
[et_target,fval] = fzero(theta_fun,t_bounds,options);

% tformat = 'YYYY-MON-DD-HR:MN:SC.####::UTC';
% Initial time
Departure_time = cspice_et2utc(et_target,'C',4);
% Compute arrival time
% Time of flight [s]
ToF = (opt_state_s(end)-opt_state_s(end-1))*TU*3600*24;
% Arrival time
et_final = et_target + ToF;
Arrival_time = cspice_et2utc(et_final,'C',4);
% PROPAGATION
% Initial conditions
z_0 = 0;
vz_0 = 0;
% Rotate to ECI
xx_rot = EM2ECI(opt_state_s(end-1), opt_state_s(1:4)', mu);
% Transform into physic unit
xx0_3D = [xx_rot(1)*DU xx_rot(2)*DU z_0 xx_rot(3)*VU xx_rot(4)*VU vz_0]';

% Celestial bodies
labels = {'Sun';
          'Mercury';
          'Venus';
          'Earth';
          'Moon';
          'Mars Barycenter';
          'Jupiter Barycenter';
          'Saturn Barycenter';
          'Uranus Barycenter';
          'Neptune Barycenter';
          'Pluto Barycenter'};
% Initialize propagation data
bodies = nbody_init(labels);
% select integration frame string (SPICE naming convention)
frame = 'J2000';
center = 'Earth';

% Propagation
options = odeset('reltol', 1e-12, 'abstol', 1e-12);
[tt_Nbody, xx_Nbody] = ode113(@(t,x) nbody_shift_rhs(t,x,bodies,frame,center), [et_target et_final], xx0_3D, options);

% Plot comparison
figure
hold on
% N-body propagator
plot3(xx_Nbody(:,1), xx_Nbody(:,2), xx_Nbody(:,3),'color','#1953D9','LineWidth',1.5);
% Simple shooting
plot3(XX(:,1)*DU,XX(:,2)*DU,zeros(length(tt),1),'LineWidth',1.5) 
plot3(0,0,0,'Marker','o','LineWidth',1,'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','b','LineStyle','none')
% Starting orbit
plot3(xx_ECI_E(1,:)*DU,xx_ECI_E(2,:)*DU,zeros(1,size(xx_ECI_E,2)),'LineStyle','-.','Color',[0.4 0.4 0.4])
% moon orbit
%plot3(rr_moon(1,:), rr_moon(2,:), rr_moon(3,:),'Color',[0.4 0.4 0.4],'LineWidth',1,'LineStyle','-.'); % from kernel
plot3(xx_Moon(1,:)*DU, xx_Moon(2,:)*DU, zeros(size(xx_Moon,2),1),'Color',[0.4 0.4 0.4],'LineWidth',1,'LineStyle','-.'); 
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
legend('N-body propagation','PBRFBP propagation','Earth','Moon orbit in the PBRFBP')
title (['@',center,':',frame])
axis equal
grid on
view(21.0592,17.3499)



% Close Kernels
cspice_kclear();

%% FUNCTIONS

% ------------------------------------------------------------------------
function XX = EM2ECI(tt, xx, mu)
%------------------------------------------------------------------------
% Rotate from Earth-Moon barycenter reference frame to ECI
%
% INPUT
% tt        [nx1] Time vector
% xx        [nx4] State vector at each time step
% mu        [1] Gravitational parameter of the attractor
%
% OUTPUT
% XX        [nx4] State vector in ECI
% ------------------------------------------------------------------------

XX(:, 1) = (xx(:, 1) + mu) .* cos(tt) - xx(:, 2) .* sin(tt);
XX(:, 2) = (xx(:, 1) + mu) .* sin(tt) + xx(:, 2) .* cos(tt);
XX(:, 3) = (xx(:,3) - xx(:,2)).* cos(tt) - (xx(:,4) + xx(:,1) + mu).* sin(tt);
XX(:, 4) = (xx(:,3) - xx(:,2)).* sin(tt) + (xx(:,4) + xx(:,1) + mu).* cos(tt);
end


%-------------------------------------------------------------------------
function DV = Delta_V (xx0_t,mu,m_s,rho,om_s, v_circ_i, v_circ_f)
%------------------------------------------------------------------------
% Objective function of a Two-Impulse Simple Shooting, computes the delta
% velocity required to perform the departure and arrival burns
%
% INPUT
% xx0_t     [(n+2)x1] State vector with final and initial time 
% mu        [1] Gravitational parameter of the main attractor
% m_s       [1] Sun scaled mass 
% rho       [1] Scaled Sun-(Earth-Moon) Batycenter distance 
% om_s      [1] Scaled ancugal velocity of the sun 
% v_circ_i  [1] Circular velocity in the departure orbit
% v_circ_f  [1] Circular velocity in the arrival orbit
%
% OUTPUT
% DV        [1] Cost of the manoeuvre
% ------------------------------------------------------------------------   

           xx_i = xx0_t(1:4) ;
           t0 = xx0_t(5);
           tf = xx0_t(6);

          [xx_f,~, ~, ~]  = propagate(t0,xx_i,tf,mu,m_s,rho,om_s);
% first delta velocity
          DV_1 = sqrt( (xx_i(3) - xx_i(2))^2 + (xx_i(4) + xx_i(1) + mu).^2  ) - v_circ_i;
% second delta velocity
          DV_2 = sqrt( (xx_f(3) - xx_f(2))^2 + (xx_f(4) + xx_f(1) + mu -1).^2  ) - v_circ_f;
% total delta velocity
          %DV = abs(DV_1) + abs(DV_2); % capire se va tenuto o tolto il valore assoluto
          DV = DV_1 + DV_2;
end

%--------------------------------------------------------------------------
function [xf,tf, xx, tt]  = propagate(t0,x0,tf,mu,m_s,rho,om_s)
%------------------------------------------------------------------------
% Orbit propagator in the PBRFBP from a starting time t0 to e final time tf
%
% INPUT
% t0        [1] Initial time
% x0        [nx1] State vector 
% tf        [1] Final time
% mu        [1] Gravitational parameter of the main attractor
% m_s       [1] Sun scaled mass 
% rho       [1] Scaled Sun-(Earth-Moon) Batycenter distance 
% om_s      [1] Scaled ancugal velocity of the sun 
%
% OUTPUT
% xf        [4x1] Final state
% tf        [1] Final time
% xx        [nx4] State eolution over time
% tt        [nx1] Time vector
% ------------------------------------------------------------------------   
   

    % Perform integration
    options = odeset('reltol', 1e-12, 'abstol', 1e-12);
    [tt, xx] = ode78(@(t,x) xyCR4BP(t,x,mu,m_s,rho,om_s), [t0 tf], x0, options);

    % Extract state vector and State Transition Matrix
    xf = xx(end,:)';
    tf = tt(end);

end

%--------------------------------------------------------------------------
function [dxdt] = xyCR4BP(t,xx, mu,ms,rho,om_s)
%------------------------------------------------------------------------
% ODE function to implement generate the equation of motion in the Planar
% Circular Restricted Four-Body Problem 
%
% INPUT
% t         [1] Time 
% xx        [4x1] State of the body ( rx, ry, vx, vz, STM(:))
% ms        [1] Sun scaled mass [-]
% rho       [1] Scaled Sun-(Earth-Moon) Batycenter distance [-]
% om_s      [1] Scaled ancugal velocity of the sun [-]
% mu        [1]  Earth–Moon mass parameter [-]
%
% OUTPUT
% dxdt      [4x1] Derivative of the state 
%------------------------------------------------------------------------
    % Extract variables
    x  = xx(1);
    y  = xx(2);
    vx = xx(3);
    vy = xx(4);

    % Potential gradient

    dUdx = x - (ms*cos(om_s*t))/rho^2 - (mu*(mu + x - 1))/((mu + x - 1)^2 + y^2)^(3/2) + ((2*mu + 2*x) * (mu - 1))/(2*((mu + x)^2 + y^2)^(3/2)) - (ms*(2*x - 2*rho*cos(om_s*t)))/(2*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(3/2));
    dUdy = y - (ms*sin(om_s*t))/rho^2 - (mu*y)/((mu + x - 1)^2 + y^2)^(3/2) - (ms*(2*y - 2*rho*sin(om_s*t)))/(2*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(3/2)) + (y*(mu - 1))/((mu + x)^2+y^2)^(3/2);

    
    % Assemble right-hand side
    dxdt = zeros(4,1);
    dxdt(1:2) = xx(3:4);
    dxdt(3)   = dUdx + 2*vy;
    dxdt(4)   = dUdy - 2*vx;
    
end

%--------------------------------------------------------------------------
function [c,ceq] = constraints(x,mu,m_s,rho,om_s,r_i,r_f)
%------------------------------------------------------------------------
% Non-Linear constraint of a Two-Impulse Simple Shooting with prescribed
% departure and arrival orbits
%
% INPUT
% x         [(4+2)x1] State vector with final and initial time 
% mu        [1] Gravitational parameter of the main attractor
% m_s       [1] Sun scaled mass 
% rho       [1] Scaled Sun-(Earth-Moon) Batycenter distance 
% om_s      [1] Scaled ancugal velocity of the sun 
% r_i       [1] Departure orbit radius
% r_f       [1] Arrival orbit radius
%
% OUTPUT
% c         [4x1] Non-linear inequality constraints
% c_eq      []    Non-linear equality constraints
% ------------------------------------------------------------------------  

% Initialization
xx_i  = x(1:4);
t0 = x(5);
tf = x(6);
% Propagation
[xx_f,~, ~, ~]  = propagate(t0,xx_i,tf,mu,m_s,rho,om_s);

% EQUALITY
% Initial: xx_i = [ x_i , y_i , vx_i, vy_i]
phi_i =  [ (xx_i(1) + mu)^2 + xx_i(2)^2 - r_i^2;
           (xx_i(1) + mu)*(xx_i(3)-xx_i(2)) + xx_i(2)*(xx_i(4)+xx_i(1) + mu)];

% Final: xx_f = [ x_f , y_f , vx_f, vy_f]
phi_f = [ (xx_f(1) + mu -1)^2 + xx_f(2)^2 - r_f^2;
          (xx_f(1) + mu -1)*(xx_f(3)- xx_f(2)) + xx_f(2)*(xx_f(4) + xx_f(1)+ mu -1)];

ceq = [phi_i;phi_f];
c = [];
end

%--------------------------------------------------------------------------
function [xf,PHIf,tf, xx, tt]  = propagate_STM(t0,x0,tf,mu,m_s,rho,om_s)
%------------------------------------------------------------------------
% Orbit and STM propagator in the PBRFBP from a starting time t0 to e final time
% tf
%
% INPUT
% t0        [1] Initial time
% x0        [4x1] State vector 
% tf        [1] Final time
% mu        [1] Gravitational parameter of the main attractor
% m_s       [1] Sun scaled mass 
% rho       [1] Scaled Sun-(Earth-Moon) Batycenter distance 
% om_s      [1] Scaled ancugal velocity of the sun 
%
% OUTPUT
% xf        [4x1] Final state
% PHIf      [4x4] STM from t0 to tf
% tf        [1] Final time
% xx        [nx4] State eolution over time
% tt        [nx1] Time vector
% ------------------------------------------------------------------------   
    % Initialize State Transition Matrix at t0
    Phi0 = eye(4);

    % Append to initial conditions the conditions for the STM
    x0Phi0 = [x0'; Phi0(:)];         
    
    % Perform integration

    % options
    options_STM = odeset('reltol', 1e-12, 'abstol', 1e-12);
    
    % ODE solver
    [tt, xx] = ode78(@(t,x) xyCR4BP_STM(t,x,mu,m_s,rho,om_s), [t0 tf], x0Phi0, options_STM);   % Call xyzCR4BP_STM  

    % Extract state vector and State Transition Matrix
    xf = xx(end,1:4)';
    PHIf = reshape(xx(end,5:end),4,4);     % from vector to STM
    tf = tt(end);

end

%--------------------------------------------------------------------------
function [dxdt] = xyCR4BP_STM(t,xx, mu,ms,rho,om_s)
%------------------------------------------------------------------------
% ODE function to implement generate the equation of motion in the Planar
% Circular Restricted Four-Body Problem and evaluate the STM
%
% INPUT
% t         [1] Time 
% xx        [20x1] State of the body ( rx, ry, vx, vz, STM(:))
% ms        [1] Sun scaled mass [-]
% rho       [1] Scaled Sun-(Earth-Moon) Batycenter distance [-]
% om_s      [1] Scaled ancugal velocity of the sun [-]
% mu        [1]  Earth–Moon mass parameter [-]
%
% OUTPUT
% dxdt      [20x1] Derivative of the state 
%------------------------------------------------------------------------

    % Extract variables
    x  = xx(1);
    y  = xx(2);
    vx = xx(3);
    vy = xx(4);
    % Put PHI in matrix form
    Phi = reshape(xx(5:end),4,4);
    % Potential gradient
    dUdx = x - (ms*cos(om_s*t))/rho^2 - (mu*(mu + x - 1))/((mu + x - 1)^2 + y^2)^(3/2) + ((2*mu + 2*x) * (mu - 1))/(2*((mu + x)^2 + y^2)^(3/2)) - (ms*(2*x - 2*rho*cos(om_s*t)))/(2*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(3/2));
    dUdy = y - (ms*sin(om_s*t))/rho^2 - (mu*y)/((mu + x - 1)^2 + y^2)^(3/2) - (ms*(2*y - 2*rho*sin(om_s*t)))/(2*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(3/2)) + (y*(mu - 1))/((mu + x)^2+y^2)^(3/2);
    % Assemble the matrix A(t)=dfdx 4x4 matrix
    dfdx =[                                                                                                                                                                                                                                                                                                                                                           0,                                                                                                                                                                                                                                                                                                                       0,  1, 0;
                                                                                                                                                                                                                                                                                                                                                        0,                                                                                                                                                                                                                                                                                                                       0,  0, 1;
(mu - 1)/((mu + x)^2 + y^2)^(3/2) - mu/((mu + x - 1)^2 + y^2)^(3/2) - ms/((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(3/2) - (3*(mu + x)^2*(mu - 1))/((mu + x)^2 + y^2)^(5/2) + (3*ms*(x - rho*cos(om_s*t))^2)/((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(5/2) + (3*mu*(2*mu + 2*x - 2)*(mu + x - 1))/(2*((mu + x - 1)^2 + y^2)^(5/2)) + 1,                                                                                (3*ms*(2*x - 2*rho*cos(om_s*t))*(2*y - 2*rho*sin(om_s*t)))/(4*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(5/2)) + (3*mu*y*(mu + x - 1))/((mu + x - 1)^2 + y^2)^(5/2) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(5/2)),  0, 2;
                                                                                                            (3*ms*(2*x - 2*rho*cos(om_s*t))*(2*y - 2*rho*sin(om_s*t)))/(4*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(5/2)) + (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(5/2)), (mu - 1)/((mu + x)^2 + y^2)^(3/2) - mu/((mu + x - 1)^2 + y^2)^(3/2) - ms/((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(3/2) - (3*y^2*(mu - 1))/((mu + x)^2 + y^2)^(5/2) + (3*ms*(y - rho*sin(om_s*t))^2)/((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(5/2) + (3*mu*y^2)/((mu + x - 1)^2 + y^2)^(5/2) + 1, -2, 0];

    % Compute the derivative of the STM
    Phidot = dfdx*Phi;
    % Assemble right-hand side
    dxdt = zeros(20,1);

    dxdt(1:2) = xx(3:4);
    dxdt(3)   = dUdx + 2*vy;
    dxdt(4)   = dUdy - 2*vx;
    dxdt(5:end) = Phidot(:);   
end

%--------------------------------------------------------------------------
function [DV,grad] = Delta_V_STM (xx0_t,mu,m_s,rho,om_s, v_circ_i, v_circ_f)
%------------------------------------------------------------------------
% Objective function of a Two-Impulse Simple Shooting evaluating the analytical derivatives,
% computes the delta velocity required to perform the departure and arrival burns
%
% INPUT
% xx0_t     [6x1] State vector with final and initial time 
% mu        [1] Gravitational parameter of the main attractor
% m_s       [1] Sun scaled mass 
% rho       [1] Scaled Sun-(Earth-Moon) Batycenter distance 
% om_s      [1] Scaled ancugal velocity of the sun 
% v_circ_i  [1] Circular velocity in the departure orbit
% v_circ_f  [1] Circular velocity in the arrival orbit
%
% OUTPUT
% DV        [1] Cost of the manoeuvre
% grad      [6x1] Derivatives of the objective function with respect to the
%                 problem variambles
% ------------------------------------------------------------------------      

           xx_i(1:4) = xx0_t(1:4) ;
           t0 = xx0_t(5);
           tf = xx0_t(6);                  
          [xx_f,PHIf,tf, ~, ~]  = propagate_STM(t0,xx_i,tf,mu,m_s,rho,om_s);
% first delta velocity
          DV_1 = sqrt( (xx_i(3) - xx_i(2))^2 + (xx_i(4) + xx_i(1) + mu).^2  ) - v_circ_i;
% second delta velocity
          DV_2 = sqrt( (xx_f(3) - xx_f(2))^2 + (xx_f(4) + xx_f(1) + mu -1).^2  ) - v_circ_f;
% total delta velocity
          DV = DV_1 + DV_2;

% GRADIENT
% Notation 
  x_i = xx_i(1);  y_i = xx_i(2); vx_i = xx_i(3); vy_i = xx_i(4); 
  x_f = xx_f(1);  y_f = xx_f(2); vx_f = xx_f(3); vy_f = xx_f(4); 
% Compute  d(DeltaV1)/d(x_i) and d(DeltaV2)/d(x_f)
D1 = [(mu + vy_i + x_i), -(vx_i - y_i), (vx_i - y_i), (mu + vy_i + x_i)]./((vx_i - y_i)^2 + (mu + vy_i + x_i)^2)^(1/2); %d(DeltaV1)/d(x_i)
D2 = [(mu + vy_f + x_f -1), -(vx_f - y_f), (vx_f - y_f), (mu + vy_f + x_f-1)]./((vx_f - y_f)^2 + (mu + vy_f + x_f-1)^2)^(1/2);  %d(DeltaV2)/d(x_f)
 
% Compute initial and final state derivative
[dxdt_i] = xyCR4BP(t0,xx_i(1:4), mu,m_s,rho,om_s);
[dxdt_f] = xyCR4BP(tf,xx_f(1:4), mu,m_s,rho,om_s);

grad = [D1 + D2*(PHIf),...
        -D2*(PHIf)*dxdt_i,...
         D2*dxdt_f           ]';
end

%--------------------------------------------------------------------------
function [c,ceq,grad,grad_ceq] = constraints_STM(x,mu,m_s,rho,om_s,r_i,r_f)
%------------------------------------------------------------------------
% Non-Linear constraint of a Two-Impulse Simple Shooting with prescribed
% departure and arrival orbits, and their derivatives with respect to the
% optimization variables
%
% INPUT
% x         [(4+2)x1] State vector with final and initial time 
% mu        [1] Gravitational parameter of the main attractor
% m_s       [1] Sun scaled mass 
% rho       [1] Scaled Sun-(Earth-Moon) Batycenter distance 
% om_s      [1] Scaled ancugal velocity of the sun 
% r_i       [1] Departure orbit radius
% r_f       [1] Arrival orbit radius
%
% OUTPUT
% c         [4x1] Non-linear inequality constraints
% c_eq      []    Non-linear equality constraints
% c         [4x1] Non-linear inequality constraints derivatives
% c_eq      []    Non-linear equality constraints derivatives
% ------------------------------------------------------------------------  


xx_i(1:4)  = x(1:4);
t0 = x(5);
tf = x(6);
[xx_f,PHIf,tf,~, ~]  = propagate_STM(t0,xx_i,tf,mu,m_s,rho,om_s);

% EQUALITY
% Initial: xx_i = [ x_i , y_i , vx_i, vy_i]
phi_i =  [ (xx_i(1) + mu)^2 + xx_i(2)^2 - r_i^2;
           (xx_i(1) + mu)*(xx_i(3)-xx_i(2)) + xx_i(2)*(xx_i(4)+xx_i(1) + mu)];
% Final: xx_f = [ x_f , y_f , vx_f, vy_f]
phi_f = [ (xx_f(1) + mu -1)^2 + xx_f(2)^2 - r_f^2;
          (xx_f(1) + mu -1)*(xx_f(3)- xx_f(2)) + xx_f(2)*(xx_f(4) + xx_f(1)+ mu -1)];
ceq = [phi_i;phi_f];
c = [];

% Gradient of the contraints
% Notation 
  x_i = xx_i(1);  y_i = xx_i(2); vx_i = xx_i(3); vy_i = xx_i(4); 
  x_f = xx_f(1);  y_f = xx_f(2); vx_f = xx_f(3); vy_f = xx_f(4); 
% Compute the derivative d(phi_i)/d(x_i) and d(phi_f)/d(x_f)
DC_i = [2*mu + 2*x_i, 2*y_i,      0, 0;      %d(phi_i)/d(x_i)
             vx_i,  vy_i, mu + x_i, y_i];
DC_f = [2*mu + 2*x_f - 2, 2*y_f,   0, 0;     %d(phi_f)/d(x_f)
            vx_f,  vy_f, mu + x_f - 1, y_f];
% Compute initial and final state time derivative
[dxdt_i] = xyCR4BP(t0,xx_i(1:4), mu,m_s,rho,om_s);
[dxdt_f] = xyCR4BP(tf,xx_f, mu,m_s,rho,om_s);

% Compute gradient of the constraints on the arrival 
g_f_xi = DC_f*PHIf;
g_f_ti = -DC_f*PHIf*dxdt_i;
g_f_tf = DC_f*dxdt_f;

grad_ceq = [DC_i(1,:), 0, 0;
            DC_i(2,:), 0, 0;
        g_f_xi(1,:), g_f_ti(1), g_f_tf(1);
        g_f_xi(2,:), g_f_ti(2), g_f_tf(2)]';
grad = [];
end

%--------------------------------------------------------------------------
function [DV,grad_V] = Delta_V_ms (xx0_t,mu, v_circ_i, v_circ_f,n_var)
%------------------------------------------------------------------------
% Objective function of a Two-Impulse Multiple Shooting with N nodes evaluating the analytical derivatives,
% computes the delta velocity required to perform the departure and
% arrival burns.
%
% INPUT
% xx0_t     [(4N+2)x1] State vector with final and initial time 
% mu        [1] Gravitational parameter of the main attractor
% v_circ_i  [1] Circular velocity in the departure orbit
% v_circ_f  [1] Circular velocity in the arrival orbit
% n_var     [1] Number of nodes N
%
% OUTPUT
% DV        [1] Cost of the manoeuvre
% grad      [(4N+2)x1] Derivatives of the objective function with respect to the
%                 problem variambles
% ------------------------------------------------------------------------    
           
           % xi, x2,....xN
           xx_i = xx0_t(1:4*n_var) ;
           xx_f = xx0_t(end-5:end-2);
% first delta velocity
          DV_1 = sqrt( (xx_i(3) - xx_i(2))^2 + (xx_i(4) + xx_i(1) + mu).^2  ) - v_circ_i;
% second delta velocity
          DV_2 = sqrt( (xx_f(3) - xx_f(2))^2 + (xx_f(4) + xx_f(1) + mu -1).^2  ) - v_circ_f;
% total delta velocity
          DV = DV_1 + DV_2;

% GRADIENT
% Notation 
  x_i = xx_i(1);  y_i = xx_i(2); vx_i = xx_i(3); vy_i = xx_i(4); 
  x_f = xx_f(1);  y_f = xx_f(2); vx_f = xx_f(3); vy_f = xx_f(4); 
% Compute  d(DeltaV1)/d(x_i) and d(DeltaV2)/d(x_f)
D1 = [(mu + vy_i + x_i), -(vx_i - y_i), (vx_i - y_i), (mu + vy_i + x_i)]./((vx_i - y_i)^2 + (mu + vy_i + x_i)^2)^(1/2); %d(DeltaV1)/d(x_i)
D2 = [(mu + vy_f + x_f -1), -(vx_f - y_f), (vx_f - y_f), (mu + vy_f + x_f-1)]./((vx_f - y_f)^2 + (mu + vy_f + x_f-1)^2)^(1/2);  %d(DeltaV2)/d(x_f)

% Gradient assembly
grad_dn = zeros(1,4*(n_var-2));
grad_V = [D1,...
        grad_dn,...
        D2,...
        0,...
        0          ]';
end

%--------------------------------------------------------------------------
function [c,ceq,d_c,d_ceq] = constraints_ms(x,mu,m_s,rho,om_s,r_i,r_f,n_var,E_RADII,M_RADII)
%------------------------------------------------------------------------
% Non-Linear constraint of a Two-Impulse Multiple Shooting with prescribed
% departure and arrival orbits, and their derivatives with respect to the
% optimization variables
%
% INPUT
% x         [(4n+2)x1] State vector with final and initial time 
% mu        [1] Gravitational parameter of the main attractor
% m_s       [1] Sun scaled mass 
% rho       [1] Scaled Sun-(Earth-Moon) Batycenter distance 
% om_s      [1] Scaled ancugal velocity of the sun 
% r_i       [1] Departure orbit radius
% r_f       [1] Arrival orbit radius
% n_var     [1] Number of nodes N
% E_RADII   [1] Earth mean radius
% M_RADII   [1] Moon mean radius
%
% OUTPUT
% c         [2nx1]       Non-linear inequality constraints
% c_eq      [4nx1]       Non-linear equality constraints
% d_c       [(4n+2)x2n]  Non-linear inequality constraints derivatives
% d_ceq     [(4n+2)x4n]  Non-linear equality constraints derivatives
% ------------------------------------------------------------------------  

% Initialization
xx_i  = x(1:4*n_var);
tf = x(end);
t0 = x(end-1);
% Time discretization
h = (tf - t0) / (n_var-1);
t_vec = t0:h:tf;
cons = NaN(4,n_var-1); 
ineq_cons = NaN(2,n_var); 

% Equality constraint gradient
grad_eq = zeros(4*n_var,4*n_var+2);
% Inequality constraint gradient
grad_ineq = zeros(2*n_var,4*n_var+2);

% Initial: xx_i = [ x_i , y_i , vx_i, vy_i]
phi_i =  [ (xx_i(1) + mu)^2 + xx_i(2)^2 - r_i^2;
           (xx_i(1) + mu)*(xx_i(3)-xx_i(2)) + xx_i(2)*(xx_i(4)+xx_i(1) + mu)];

% n-constraints
for j = 1 : length(t_vec)-1
    ind = (1:4)+ 4*(j-1);
    x_j = x(ind);
     % INEQUALITY CONSTRAINT
    ineq_cons(:,j) = [ E_RADII(1)^2 - (x_j(1)+ mu)^2 - x_j(2)^2
                       M_RADII(1)^2 - (x_j(1)+ mu-1)^2 - x_j(2)^2 ];
     % GRAD INEQUALITY
    it = (1:2)+ 2*(j-1);
    S = [ -2*(x_j(1) + mu), -2*x_j(2), 0, 0;
           -2*(x_j(1) + mu-1), -2*x_j(2), 0, 0];
    grad_ineq(it,ind) = S;

    % EQUALITY CONSTRAINT
    [x_j,PHIf,~, ~, ~]  = propagate_STM(t_vec(j),x_j',t_vec(j+1),mu,m_s,rho,om_s);
    cons(:,j) = x_j- x(ind+4);
    % Compile diagonal and upper diagonal terms from 1 to N-1
    grad_eq(ind,ind) = PHIf;
    grad_eq(ind,ind+4) = -eye(4);

    % Compute f(x_j, t_j)
    [dxdt_j] = xyCR4BP(t_vec(j),x(ind), mu,m_s,rho,om_s);
    % Compute f( phi(x_j,t_j;t_j+1), t_j+1)
    [dxdt_j_plus] = xyCR4BP(t_vec(j+1),x_j, mu,m_s,rho,om_s);
    % assamble to find Qj-1, Qj-N
    Q_j_1 = (j-n_var)/(n_var-1)*PHIf*dxdt_j + (n_var-j-1)/(n_var-1)*dxdt_j_plus;
    Q_j_N = (1-j)/(n_var-1)*PHIf*dxdt_j     +  j/(n_var-1)*dxdt_j_plus         ;
    grad_eq(ind,end-1:end) = [Q_j_1 , Q_j_N];
end
% x final
xx_f = x(end-5:end-2);

% GRADIENT of final and intial constraints
   grad_eq(end-3:end-2,1:4) = [2*mu + 2*xx_i(1), 2*xx_i(2),      0, 0;      %d(phi_i)/d(x_i)
                               xx_i(3),  xx_i(4), mu + xx_i(1), xx_i(2)];

   grad_eq(end-1:end,end-5:end-2)= [2*mu + 2*xx_f(1)-2, 2*xx_f(2),      0, 0;      %d(phi_f)/d(x_n)
                               xx_f(3),  xx_f(4), mu + xx_f(1)-1, xx_f(2)];

   % Gradient of final ineq constraints
   S = [ -2*(xx_f(1) + mu), -2*xx_f(2), 0, 0;
           -2*(xx_f(1) + mu-1), -2*xx_f(2), 0, 0]; 
    grad_ineq(end-1:end,end-5:end-2) = S;
% Last inequality constraint
ineq_cons(:,end) = [ E_RADII(1)^2 - (xx_f(1)+ mu)^2 - xx_f(2)^2
                   M_RADII(1)^2 - (xx_f(1)+ mu-1)^2 - xx_f(2)^2  ];
% Final equality constraint: xx_f = [ x_f , y_f , vx_f, vy_f]
phi_f = [ (xx_f(1) + mu -1)^2 + xx_f(2)^2 - r_f^2;
          (xx_f(1) + mu -1)*(xx_f(3)- xx_f(2)) + xx_f(2)*(xx_f(4) + xx_f(1)+ mu -1)];

% EQUALITY CONSTRAINT
ceq = [cons(:);phi_i; phi_f];
% INEQUALITY CONSTRAINT
c = ineq_cons(:);
% EQUALITY CONSTRAINT GRADIENT
d_ceq = grad_eq';
% INEQUALITY CONSTRAINT GRADIENT
d_c = grad_ineq' ;
end

%-------------------------------------------------------------------------- 
function theta_res = Thetafinder(t,theta_i)
%-------------------------------------------------------------------------
% Objective function for the epoch finding, compare the target
% Earth-Moon-Sun angle with the one computed from kernel at each time,
% the function compute the residual.
%
% INPUT
% t         [1] Time instant
% theta_i   [1] Target EMS angle
%
% OUTPUT
% theta_res [1] Residual 
%------------------------------------------------------------------------
    
    % Moon position in J2000
    Moon_pos = cspice_spkpos('MOON',t,'J2000','NONE','EARTH');
    % SUN position in J2000
    Sun_pos = cspice_spkpos('SUN',t,'J2000','NONE','EARTH');

    theta = atan2(Sun_pos(2,:),Sun_pos(1,:)) - atan2(Moon_pos(2,:),Moon_pos(1,:)) ;
    theta_res = wrapTo2Pi(theta) - theta_i;
end



%--------------------------------------------------------------------------
function [bodies] = nbody_init(labels)
%NBODY_INIT Initialize planetary data for n-body propagation
%   Given a set of labels of planets and/or barycentres, returns a
%   cell array populated with structures containing the body label and the
%   associated gravitational constant.
%
%
% Author
%   Name: ALESSANDRO 
%   Surname: MORSELLI
%   Research group: DART
%   Department: DAER
%   University: Politecnico di Milano 
%   Creation: 26/09/2021
%   Contact: alessandro.morselli@polimi.it
%   Copyright: (c) 2021 A. Morselli, Politecnico di Milano. 
%                  All rights reserved.
%
%
% Notes:
%   This material was prepared to support the course 'Satellite Guidance
%   and Navigation', AY 2021/2022.
%
%
% Inputs:
%   labels : [1,n] cell-array with object labels
%
% Outputs:
%   bodies : [1,n] cell-array with struct elements containing the following
%                  fields
%                  |
%                  |--bodies{i}.name -> body label
%                  |--bodies{i}.GM   -> gravitational constant [km**3/s**2]
%
%
% Prerequisites:
%   - MICE (Matlab SPICE)
%   - Populated kernel pool (PCK kernels)
%

% Initialize output
bodies = cell(size(labels));

% Loop over labels
for i = 1:length(labels)
    % Store body label
    bodies{i}.name = labels{i};
    % Store body gravitational constant
    bodies{i}.GM   = cspice_bodvrd(labels{i}, 'GM', 1);
end

end

%--------------------------------------------------------------------------
function [dxdt] = nbody_shift_rhs(t, x, bodies, frame, center)
% ---------------------------------------------------------------------
%   Evaluates the right-hand-side of a newtonian N-body propagator.
%   
% INPUT
% t       [1,1] ephemeris time (ET SPICE), seconds past J2000 (TDB)
% x       [6,1] cartesian state vector wrt desired object
% bodies  [1,n] cell-array created with function nbody_init
% center  [str] string with the choosen center
%
% OUTPUT
% dxdt    [6,1] RHS, newtonian gravitational acceleration only
%-------------------------------------------------------------------------

if not( strcmpi(frame, 'ECLIPJ2000') || strcmpi(frame, 'J2000') )
    msg = 'Invalid integration reference frame, select either J2000 or ECLIPJ2000';
    error(msg);
end
% Initialize right-hand-side
dxdt = zeros(6,1);
% Position derivative is object's velocity
dxdt(1:3) = x(4:6);
% Extract the object position from state x
rr_center_obj = x(1:3);
% Extract GM of central body GM0
GM0 = cspice_bodvrd(center, 'GM', 1);
 % Compute square distance and distance
dist_center_2 = dot(rr_center_obj, rr_center_obj);
dist_center = sqrt(dist_center_2);
    % Compute the gravitational acceleration using Newton's law
aa_grav_GM0 =  - GM0* rr_center_obj /(dist_center*dist_center_2);
% Loop over all bodies (except GM0)
for i = 1:length(bodies)
    if ~strcmpi(bodies{i}.name,center)

    % rho
    rv_center_body = cspice_spkezr(bodies{i}.name, t, frame, 'NONE', center);
    rho = rv_center_body(1:3);
    % Extract object position wrt. i-th celestial body (d)
    rr_body_obj = rr_center_obj - rv_center_body(1:3);
    dist2 = dot(rr_body_obj, rr_body_obj);
    dist = sqrt(dist2);
    % Compute (q, f):
    q = dot( rr_center_obj, (rr_center_obj-2*rho) ) / dot(rho,rho);
    f = q*(3+3*q+q^2) / (1+(1+q)^(3/2));
    % Compute their contribution to acceleration:
    acc = - bodies{i}.GM / (dist2*dist)*(rr_center_obj + rho*f);
    % Sum up acceleration to right-hand-side
    dxdt(4:6) = dxdt(4:6) + acc;
    end
end
 dxdt(4:6) = dxdt(4:6) + aa_grav_GM0;
end

%--------------------------------------------------------------------------
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

