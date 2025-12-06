% Spacecraft Guidance and Navigation (2024/2025)
% Assignment # 2, Exercise 1
% Author: Davide Bellini

clearvars; close all; clc;
cspice_kclear()

plot_set(); % Plot settings

% Upload kernels
% Load kernels
cspice_furnsh('assignment02.tm');

%% DATA
% DATA
% Iinitial state 
r_i = [-0.011965533749906, -0.017025663128129];
v_i = [10.718855256727338,  0.116502348513671];

% Initial time 
t_i = 1.282800225339865;
t_f = 9.595124551366348;

% Covariance
P0 = [ 1.041e-15  6.026e-17  5.647e-16  4.577e-15;
      6.026e-17  4.287e-18  4.312e-17  1.855e-16;
      5.647e-16  4.312e-17  4.432e-16  1.455e-15;
      4.577e-15  1.855e-16  1.455e-15  2.822e-14  ];

% Data for UT
alpha = 1;
beta = 2; 
k = 0;

% DATA from ∗F. Topputo, “On optimal two-impulse Earth–Moon transfers in a four-body model”, Celestial Mechanics and Dynamical Astronomy, Vol. 117, pp. 279–313, 2013
m_sun = 3.28900541*1e5;   % [-] Sun scaled mass
rho = 3.88811143*1e2;     % [-] Scaled Sun-(Earth-Moon) Batycenter distance
om_s = -9.25195985*1e-1;  % [-] Scaled ancugal velocity of the sun
 mu =  1.21506683*1e-2;   % [-]  Earth–Moon mass parameter

%% ------------------ 1.1 Mean and Covariance propagation ----------------

% Time vector
n_points = 5;
timePoints = linspace(t_i,t_f,n_points);

% Initialization;
xx_0 = [r_i, v_i];        % Initial state
LinCov = struct('mean', {}, 'covariance', {});
UT = struct('mean', {}, 'covariance', {});
MC = struct('mean', {}, 'covariance', {});

LinCov(1).mean = xx_0';
LinCov(1).covariance = P0;
UT(1).mean = xx_0';
UT(1).covariance = P0;
MC(1).mean = xx_0';
MC(1).covariance = P0;

% Pre-Processing for UT
    n = length(xx_0);
    lambda = alpha^2*(n + k) - n;
    sigmaPoints = zeros(n,1+2*n);
    sqrt_Mat = chol((n+lambda)*P0,"lower");
    propSigmaPoints =  zeros(n,1+2*n);

% Sigma Points generation
    sigmaPoints(:,1) =  xx_0';
    sigmaPoints(:,2:n+1) = repmat(xx_0', 1, n) + sqrt_Mat;  
    sigmaPoints(:,n+2:2*n+1) = repmat(xx_0', 1, n) - sqrt_Mat; 

    % Weights
    W_m0 = lambda / (lambda + n);
    W_c0 = lambda / (lambda + n) + (1 - alpha^2 + beta);
    Wcm = 1/ (2*(n+lambda));
    W_covariance = [W_c0; Wcm*ones(2*n,1)]';
    W_mean = [W_m0; Wcm*ones(2*n,1)]';

%--------------------- Linearized Approach (LinCov)----------------------
for i = 1:n_points-1
    
    t_end = timePoints(i+1);
    % Propagation
    [xf,PHIf,tf, xx, ~]  = propagate_STM(t_i,xx_0,t_end,mu,m_sun,rho,om_s);
    
    % Mean and Covariance
    LinCov(i+1).mean = xf;
    LinCov(i+1).covariance = PHIf*P0*PHIf';

end

%---------------------Unscented Transform (UT)-------------------------
for i = 1: n_points-1
    % Time interval
    t_end = timePoints(i+1);

    for j = 1: size(sigmaPoints,2)
        % Propagation of sigma points
        [yyf,~,~, ~, ~]  = propagate_STM(t_i, sigmaPoints(:,j)',t_end,mu,m_sun,rho,om_s);
        propSigmaPoints(:,j) = yyf;
    end

     % Mean and Covariance
     UT(i+1).mean = sum(W_mean.*propSigmaPoints,2);
     UT(i+1).covariance = (propSigmaPoints - repmat(UT(i+1).mean, 1, 9))*diag(W_covariance)*(propSigmaPoints - repmat(UT(i+1).mean, 1, 9))';

end

% Level of confidence
num_sigma = 3;

% Creating the circle for plot
theta = linspace(0, 2*pi, 1000);
circle = [cos(theta); sin(theta)];

% Eigenvalues and eigenvector of the covariance matrix LC
[V_LC, D_LC] = eig(LinCov(end).covariance(1:2,1:2));
% Create the scaled ellipse
ellipse_LC = V_LC * sqrt(D_LC) * (num_sigma * circle);
% Ellipse traslation in the mean
ellipse_LC = ellipse_LC + LinCov(end).mean(1:2);

% Eigenvalues and eigenvector of the covariance matrix UT
[V_UT, D_UT] = eig(UT(end).covariance(1:2,1:2));
% Create the scaled ellipse
ellipse_UT = V_UT * sqrt(D_UT) * (num_sigma * circle);
% Ellipse traslation in the mean
ellipse_UT = ellipse_UT + UT(end).mean(1:2);


% Plot
figure;
hold on;
plot(ellipse_LC(1, :), ellipse_LC(2, :), 'Color',[0.82,0.41,0.12], 'LineWidth', 1.5,'LineStyle','-.');  
plot(LinCov(end).mean(1), LinCov(end).mean(2),'o','MarkerFaceColor', [0.82,0.41,0.12],'MarkerSize',8); 
plot(ellipse_UT(1, :), ellipse_UT(2, :),'Color',[0.27,0.35,0.47], 'LineWidth',1.5,'LineStyle','-');  
plot(UT(end).mean(1), UT(end).mean(2),'s','MarkerFaceColor',[0.27,0.35,0.47],'MarkerSize',8,'MarkerEdgeColor',[0.27,0.35,0.47]); 
title('Mean Position and Confidence Ellipse (3$\sigma$) at final time');
xlabel('$x$ [-]');
ylabel('$y$ [-]');
legend('Covariance LinCov','Mean LinCov','Covariance UT','Mean UT')
grid on;
axis equal;


%% ------------------ 1.2 Montecarlo Simulation -----------------------

% Initialization
pop = 1000;
samples = zeros(n,pop);

% MC samples generation
R = mvnrnd(xx_0, P0, pop);
for i = 1:n_points-1
    
    t_end = timePoints(i+1);
    % samples propagation
    for j = 1:size(R,1)
    xx = R(j,:);
    [xf,PHIf,tf, xx, ~]  = propagate_STM(t_i,xx,t_end,mu,m_sun,rho,om_s);
    samples(:,j) = xf; 
    end

    % Sample mean and covariacne
    samp_mean = 1/pop * sum(samples,2);
    samp_var  = 1/ (pop-1) * (samples - repmat(samp_mean, 1, pop))*(samples - repmat(samp_mean, 1, pop))';

    MC(i+1).mean = samp_mean;
    MC(i+1).covariance = samp_var;

end

% Eigenvalues and eigenvector of the covariance matrix montecarlo
[V_MC, D_MC] = eig(MC(end).covariance(1:2,1:2));
% Create the scaled ellipse
ellipse_MC = V_MC * sqrt(D_MC) * (num_sigma * circle);
% Ellipse traslation in the mean
ellipse_MC = ellipse_MC + MC(end).mean(1:2);


% Plot MC
figure;
hold on;
plot(ellipse_LC(1, :), ellipse_LC(2, :), 'Color',[0.82,0.41,0.12], 'LineWidth', 1.5,'LineStyle','-.'); 
plot(LinCov(end).mean(1), LinCov(end).mean(2),'o','MarkerFaceColor', [0.82,0.41,0.12],'MarkerSize',8); 

plot(ellipse_UT(1, :), ellipse_UT(2, :),'Color',[0.27,0.35,0.47], 'LineWidth',1.5,'LineStyle','-');  
plot(UT(end).mean(1), UT(end).mean(2),'s','MarkerFaceColor',[0.27,0.35,0.47],'MarkerSize',8,'MarkerEdgeColor',[0.27,0.35,0.47]); 

plot(ellipse_MC(1, :), ellipse_MC(2, :),'Color',[0.3, 0.85, 0.4],'LineStyle','--', 'LineWidth',1.5);  
plot(samp_mean(1), samp_mean(2),'d','MarkerFaceColor',[0.3, 0.85, 0.4],'MarkerSize',8,'MarkerEdgeColor',[0.3, 0.85, 0.4]);
plot(samples(1,:),samples(2,:),'.','Color','k','MarkerSize',5)

title('Mean Position and Confidence Ellipse (3$\sigma$) at final time');
xlabel('$x$ [-]');
ylabel('$y$ [-]');
legend('Covariance LinCov','Mean LinCov','Covariance UT','Mean UT','Covariance MC','Mean MC','Propagated MC samples ')
grid on;
axis equal;

% Max eigenvalues
% Initialization
lambda_LinCov = zeros(2,n_points);
lambda_UT = zeros(2,n_points);
lambda_MC = zeros(2,n_points);

% Linear Covariance
for i = 1:n_points
   
    % Maximum eigenvalue: Position
    Pr = LinCov(i).covariance(1:2,1:2);
    tracePr = trace(Pr);         
    detPr = det(Pr);             
    lambda_max_r = (tracePr + sqrt(tracePr^2 - 4 * detPr)) / 2;
     % Maximum eigenvalue: Velocity
    Pv = LinCov(i).covariance(3:4,3:4);
    tracePv = trace(Pv);         
    detPv = det(Pv);             
    lambda_max_v = (tracePv + sqrt(tracePv^2 - 4 * detPv)) / 2;

    lambda_LinCov(1:2,i) = [lambda_max_r;lambda_max_v];
end

% Unscented Transform
for i = 1:n_points
   
    % Maximum eigenvalue: Position
    Pr = UT(i).covariance(1:2,1:2);
    tracePr = trace(Pr);         
    detPr = det(Pr);             
    lambda_max_r = (tracePr + sqrt(tracePr^2 - 4 * detPr)) / 2;
     % Maximum eigenvalue: Velocity
    Pv = UT(i).covariance(3:4,3:4);
    tracePv = trace(Pv);         
    detPv = det(Pv);             
    lambda_max_v = (tracePv + sqrt(tracePv^2 - 4 * detPv)) / 2;

    lambda_UT(1:2,i) = [lambda_max_r;lambda_max_v];
end

% Linear Covariance
for i = 1:n_points
   
    % Maximum eigenvalue: Position
    Pr = MC(i).covariance(1:2,1:2);
    tracePr = trace(Pr);         
    detPr = det(Pr);             
    lambda_max_r = (tracePr + sqrt(tracePr^2 - 4 * detPr)) / 2;
     % Maximum eigenvalue: Velocity
    Pv = MC(i).covariance(3:4,3:4);
    tracePv = trace(Pv);         
    detPv = det(Pv);             
    lambda_max_v = (tracePv + sqrt(tracePv^2 - 4 * detPv)) / 2;

    lambda_MC(1:2,i) = [lambda_max_r;lambda_max_v];
end

sigma_LinCov = num_sigma*sqrt(lambda_LinCov);
sigma_UT = num_sigma*sqrt(lambda_UT);
sigma_MC = num_sigma*sqrt(lambda_MC);


% Plot 3 sigma max
figure;

hold on;
plot(timePoints, sigma_LinCov(1, :), 'Color',[0.82,0.41,0.12], 'LineWidth', 1, 'Marker', 's', 'MarkerSize', 6, 'LineStyle', '-'); 
plot(timePoints, sigma_UT(1, :),'Color',[0.27,0.35,0.47], 'LineWidth',1, 'Marker', 'd', 'MarkerSize', 6, 'LineStyle', '--');  
plot(timePoints, sigma_MC(1, :),'Color',[0.3, 0.85, 0.4],'LineWidth',1, 'Marker', 'o', 'MarkerSize', 6, 'LineStyle', '-.');  
title('3$\sigma_{max}$: Position');
xlabel('$t$  [-]');
ylabel('$3 \ \sigma_{max}$ [-]');
legend('LinCov','UT','MC','Location','northwest')
grid on;
ylim([-1 20]*1e-4)

figure;
plot(timePoints, sigma_LinCov(2, :), 'Color',[0.82,0.41,0.12], 'LineWidth', 1, 'Marker', 's', 'MarkerSize', 6, 'LineStyle', '-');
hold on;
plot(timePoints, sigma_UT(2, :), 'Color',[0.27,0.35,0.47], 'LineWidth', 1, 'Marker', 'd', 'MarkerSize', 6, 'LineStyle', '--');
plot(timePoints, sigma_MC(2, :), 'Color',[0.3, 0.85, 0.4], 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 6, 'LineStyle', '-.');
title('3$\sigma_{max}$: Velocity');
xlabel('$t$ [-]');
ylabel('$3 \ \sigma_{max}$ [-]');
legend('LinCov','UT','MC','Location','northwest')
grid on;
ylim([-0.01 0.5])

% qqplot
figure;
subplot(2,2,1)
plotSet = qqplot(samples(1,:));
plotSet(1).Marker = 'x';  
plotSet(1).MarkerSize = 5;  
plotSet(1).MarkerEdgeColor = [0.35,0.35,0.47] ;
plotSet(1).LineWidth = 0.8;
plotSet(2).LineWidth = 1.7;
plotSet(3).LineWidth = 1.7;
grid on
title('Position: $x$ direction')

subplot(2,2,2)
plotSet = qqplot(samples(2,:));
plotSet(1).Marker = 'x';  
plotSet(1).MarkerSize = 5;  
plotSet(1).MarkerEdgeColor = [0.35,0.35,0.47] ;
plotSet(1).LineWidth = 0.8;
plotSet(2).LineWidth = 1.7;
plotSet(3).LineWidth = 1.7;
grid on
title('Position: $y$ direction')

subplot(2,2,3)
plotSet = qqplot(samples(3,:));
plotSet(1).Marker = 'x';  
plotSet(1).MarkerSize = 5;  
plotSet(1).MarkerEdgeColor = [0.35,0.35,0.47] ;
plotSet(1).LineWidth = 0.8;
plotSet(2).LineWidth = 1.7;
plotSet(3).LineWidth = 1.7; 
grid on
title('Velocity: $x$ direction')

subplot(2,2,4)
plotSet = qqplot(samples(4,:));
plotSet(1).Marker = 'x';  
plotSet(1).MarkerSize = 5;  
plotSet(1).MarkerEdgeColor = [0.35,0.35,0.47] ;
plotSet(1).LineWidth = 0.8;
plotSet(2).LineWidth = 1.7;
plotSet(3).LineWidth = 1.7;  
grid on
title('Velocity: $y$ direction')


%% FUNCTIONS

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

function [xf,PHIf,tf, xx, tt]  = propagate_STM(t0,x0,tf,mu,m_s,rho,om_s)
% -------------------------------------------------------------------
% Propagate an orbit from an initial time t0 to a final time tf, given the
% initial state conditions x0 and the mu constant
%
% INPUT
% t0        [1] Initial Time
% tf        [1] Final Time
% x0        [4x1] Initial state
% ms        [1] Sun scaled mass [-]
% rho       [1] Scaled Sun-(Earth-Moon) Batycenter distance [-]
% om_s      [1] Scaled ancugal velocity of the sun [-]
% mu        [1]  Earth–Moon mass parameter [-]
%
% OUTPUT
% xf        [4x1] Final state
% PHI       [4x4] STM from t0 to tf
% tf        [1] Final time
% xx        [nx4] State eolution over time
% tt        [nx1] Time vector
%---------------------------------------------------------------------
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