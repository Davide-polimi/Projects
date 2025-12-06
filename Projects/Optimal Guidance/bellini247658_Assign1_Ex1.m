% Spacecraft Guidance and Navigation (2024/2025)
% Assignment # 1, Exercise 1
% Author: Davide Bellini

clearvars, close all; clc;
% Plot settings
plot_set();

%% ----------------------------- 1.1 ---------------------------------------
% Coordinates of the five Lagrange
% mu
mu = 0.012150; 

% Potential function
dU = @(x,y,z) [ x - (1-mu).*(mu+x)./ ( sqrt( (x+mu).^2 + y.^2 + z.^2)).^3 + mu.*(1-mu-x)./ ( sqrt( (x+mu-1).^2 + y.^2 + z.^2)).^3 
                y - (1-mu).*y./ ( sqrt( (x+mu).^2 + y.^2 + z.^2)).^3 - mu.*y./ ( sqrt( (x+mu-1).^2 + y.^2 + z.^2)).^3
                 - (1-mu).*z./ ( sqrt( (x+mu).^2 + y.^2 + z.^2)).^3 - mu.*z./ ( sqrt( (x+mu-1).^2 + y.^2 + z.^2)).^3];

% Lagrangian point first guesses (from theoretical consideration on report)

% Collinear Points
L1_0 = [0.9,0,0];     % Between Earth and Moon
L2_0 = [1.2,0,0];     % After the moon
L3_0 = [-0.2,0,0];    % Opposite to L2

% Triangular Points
L4_0 = [1/2, sqrt(3)/2,0];
L5_0 = [1/2, -sqrt(3)/2,0];

% Solver option
 options = optimoptions('fsolve', 'Display', 'iter','FunctionTolerance',1e-12,'OptimalityTolerance',1e-12); 

% Collinear libration points
L1 = fsolve(@(x) dU(x(1), x(2), x(3)), L1_0, options);
L2 = fsolve(@(x) dU(x(1), x(2), x(3)), L2_0, options);
L3 = fsolve(@(x) dU(x(1), x(2), x(3)), L3_0, options);

% Triangular libration points
L4 = fsolve(@(x) dU(x(1), x(2), x(3)), L4_0, options);
L5 = fsolve(@(x) dU(x(1), x(2), x(3)), L5_0, options);

% Earth position
P1 = [-mu; 0];
% Moon position
P2 = [1-mu, 0];

% PLOT LIBRATION POINTS
% text offset for the plot
offset = 0.08;

figure;
hold on
% Plot circumferences for triangular points
theta = linspace(0, 2*pi, 100);  r = 1;
xEarth = -mu + r * cos(theta);  yEarth = r * sin(theta); 
xMoon = (1 - mu) + r * cos(theta);  yMoon = r * sin(theta); 
plot(xEarth, yEarth, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.6, 'HandleVisibility', 'off')  
plot(xMoon, yMoon, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.6, 'HandleVisibility', 'off')  
% Plot function for collinear points:
xvec = linspace(-1.7, 2.1, 1000);
dUdx = xvec - (1-mu)*(xvec+mu)./abs(xvec+mu).^3 - mu*(xvec+mu-1)./abs(xvec+mu-1).^3;
dUdx(dUdx > 15 | dUdx < -15) = NaN;         % to avoid discontinuities
plot(xvec, dUdx, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.6, 'HandleVisibility', 'off')
% Plot x-axis, y-axis
line([-1.8,1.8], [0 0], 'Color', 'k', 'LineWidth',0.5);  % x
line([0 0],[-1.8,1.8], 'Color', 'k', 'LineWidth', 0.5);  % y
% Draw the arrows at the end
quiver(-1.8, 0, 4, 0, 'Color', 'k', 'LineWidth', 0.5, 'MaxHeadSize', 0.2);
quiver(0, -1.8, 0, 4, 'Color', 'k', 'LineWidth', 0.5, 'MaxHeadSize', 0.2);
text(1.78, -0.15, 'x', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
text(-0.15, 1.78, 'y', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

% Plot Earth and Moon
plot(P1(1),P1(2),'Marker','o','LineWidth',1,'MarkerSize',15,'MarkerEdgeColor','k','MarkerFaceColor','b')
plot(P2(1),P2(2),'Marker','o','LineWidth',1,'Color',[0.7, 0.7, 0.7],'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',[0.7, 0.7, 0.7])
% Plot Libration points 
plot(-mu+L1(1),L1(2),'Marker','o','LineWidth',1,'Color','r','MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor','r')
plot(-mu+L2(1),L2(2),'Marker','o','LineWidth',1,'Color','r','MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor','r')
plot(-mu+L3(1),L3(2),'Marker','o','LineWidth',1,'Color','r','MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor','r')
plot(L4(1),L4(2),'Marker','o','LineWidth',1,'Color','r','MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor','r')
plot(L5(1),L5(2),'Marker','o','LineWidth',1,'Color','r','MarkerSize',7,'MarkerEdgeColor','k','MarkerFaceColor','r')
% Plot text on the figure
text(P1(1), P1(2)-offset, 'Earth','HorizontalAlignment','center', 'VerticalAlignment', 'cap','FontWeight','bold','FontSize',13);
text(P2(1), P2(2)-offset, 'Moon','HorizontalAlignment','center', 'VerticalAlignment', 'top','FontWeight','bold','FontSize',13);
text(L1(1), offset, 'L1','HorizontalAlignment','center', 'VerticalAlignment', 'bottom','FontWeight','bold','FontSize',13);
text(L2(1), offset, 'L2','HorizontalAlignment','center', 'VerticalAlignment', 'bottom','FontWeight','bold','FontSize',13);
text(L3(1), offset, 'L3','HorizontalAlignment','center', 'VerticalAlignment', 'bottom','FontWeight','bold','FontSize',13);
text(L4(1), L4(2) + offset, 'L4','HorizontalAlignment','center', 'VerticalAlignment', 'bottom','FontWeight','bold','FontSize',13);
text(L5(1), L5(2) + offset, 'L5','HorizontalAlignment','center', 'VerticalAlignment', 'bottom','FontWeight','bold','FontSize',13);
% Plot options
xlim([-2,2])
ylim([-2,2])
grid on
axis equal
xlabel('x [DU]')
ylabel('y [DU]')
title('Lagrangian points in the Earth-Moon system','FontSize',15)

% JACOBI constant

% Potential definition
OM = @(x,y,z) 1/2*(x.^2 + y.^2) + (1-mu)./sqrt( (x+mu).^2 +y.^2 +z.^2) + mu./sqrt( (x+mu-1).^2 +y.^2 +z.^2) +1/2*mu*(1-mu);
% Compute Jacobi Constant C = 2*OM - v^2 for each point (velocity is zero)
C1 = 2* OM(L1(1),L1(2),L1(3));
C2 = 2* OM(L2(1),L2(2),L2(3));
C3 = 2* OM(L3(1),L3(2),L3(3));
C4 = 2* OM(L4(1),L4(2),L4(3));
C5 = 2* OM(L5(1),L5(2),L5(3));

%% ----------------------------- 1.2 -------------------------
% Halo orbits
% Initial state data
x0  = 1.068792441776;
y0  = 0;
z0  = 0.071093328515;
vx0 = 0;
vy0 = 0.319422926485;
vz0 = 0;
xx0 = [ x0; y0; z0; vx0; vy0; vz0];

tf = 5.0;   % Propagation time set to a sufficient high value to ensure the y-axis crossing
[~,~,~,xx_F]  = propagate3D(0,xx0,tf,mu,false);  % set flag  to false to continue the propagation up to final time
[~,~,~,xx_B]  = propagate3D(0,xx0,-tf,mu,false); % set flag  to false to continue the propagation up to final time

% Plot of the propagated initial conditions
figure
hold on
grid on
plot3(xx_F(:,1),xx_F(:,2),xx_F(:,3),'LineWidth',2)
plot3(xx_B(:,1),xx_B(:,2),xx_B(:,3),'LineWidth',2)
% Moon
plot3(P2(1), P2(2), 0, 'Marker', 'o', 'LineWidth', 1, 'Color', 'none', ...
    'MarkerSize', 12, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.7, 0.7, 0.7], 'HandleVisibility', 'off');
text(P2(1), P2(2)-offset/4, 0, 'Moon', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
    'FontWeight', 'bold', 'FontSize', 13);
% L2
plot3(-mu+L2(1), L2(2), 0, 'Marker', 'o', 'LineWidth', 1, 'Color', 'none', ...
    'MarkerSize', 7, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
text(L2(1), L2(2), 0, 'L2', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'FontWeight', 'bold', 'FontSize', 13);
xlabel('x [DU]')
ylabel('y [DU]')
zlabel('z [DU]')
axis equal
legend(["Forward","Backward"],'Location','best')
view(14.5342,9.7729)

% Correction
% Errors inzialization
err_y = 1;    
err_vxf = 1;    
err_vzf = 1;
err_C = 1;
% Initialization
Nmax    = 100;    % Set a maximum number of iterations to avoid infinite looping
iter    = 0;      % Initialize the iteration counter
tol     = 1e-12;  % Set the desired tolerance
% First guesses
xx0_new = xx0;
MAT = zeros(4,4);
correction = [xx0_new(1), xx0_new(3), xx0_new(5), tf]';
% Target JACOBI CONSTANT 
C_ref = 3.09;
% Iterative corrections
while (abs(err_vxf) > tol || abs(err_vzf) > tol || abs(err_C) > tol || abs(err_y) > tol) && iter < Nmax

    % Evaluate the Jacobi constant and its gradient wrt new initial conditions
    [C_iter, dC_dx0] = JacobiConstant(xx0_new,mu);

    % Perform propagation up to event time
    [xf,PHI,te]  = propagate3D(0,xx0_new, tf , mu);
   
    % Evaluate the flow derivative wrt the final time EOM (dphi/dte)
    [dxdt] = xyzCR3BP_STM(0,[xf;PHI(:)], mu);   % note PHI doesn't affect the state components of the output

    % Assembling the transition matrix MAT  with elements of the previous results
    MAT(1:3,1:3) = [PHI(2,1), PHI(2,3) PHI(2,5);        % Add the STM components
                    PHI(4,1), PHI(4,3) PHI(4,5);
                    PHI(6,1), PHI(6,3) PHI(6,5)];

    MAT(4,:) =     [dC_dx0(1),dC_dx0(3) dC_dx0(5), 0];  % Add the Jacobin Constant gradient components 

    MAT(:,4) =     [dxdt(2),dxdt(4),dxdt(6),0];         % Add the time derivative components

    % Compute the deviation in the final state
      err_y =   xf(2);
      err_vxf = xf(4);
      err_vzf = xf(6);
      err_C =  C_iter -  C_ref ;
      errors = [err_y, err_vxf, err_vzf, err_C]';  % error vector     
    % Compute the correcction 
      correction = correction - MAT \ errors; 
    % Update the initial conditions
      xx0_new(1) = correction(1);
      xx0_new(3) = correction(2);
      xx0_new(5) = correction(3);
              tf = correction(4);  
    % Update iteration counter
    iter = iter+1;
end

% Propagation of the orbit
% Forward
[~,~,te, xx_F, tt_F] = propagate3D(0,xx0_new, te,mu,false);
% Backward
[~,~,~, xx_B, tt_B] = propagate3D(0,xx0_new, -te,mu,false);

% PLOT of the corrected orbit
figure 
hold on
grid on
% Halo orbit
plot3(xx_F(:,1),xx_F(:,2),xx_F(:,3),'LineWidth',2)
plot3(xx_B(:,1),xx_B(:,2),xx_B(:,3),'LineWidth',2)
% Moon
plot3(P2(1), P2(2), 0, 'Marker', 'o', 'LineWidth', 1, 'Color', 'none', ...
    'MarkerSize', 12, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.7, 0.7, 0.7], 'HandleVisibility', 'off');
text(P2(1), P2(2)-offset/4, 0, 'Moon', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
    'FontWeight', 'bold', 'FontSize', 13);
% L2
plot3(-mu+L2(1), L2(2), 0, 'Marker', 'o', 'LineWidth', 1, 'Color', 'none', ...
    'MarkerSize', 7, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
text(L2(1), L2(2), 0, 'L2', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'FontWeight', 'bold', 'FontSize', 13);
legend('Forward','Backward','Location','best')
xlabel('x [DU]')
ylabel('y [DU]')
zlabel('z [DU]')
axis equal
title('Halo orbit ( C = 3.09) ','FontSize',15)
view(14.5342,9.7729)

%% -------------------------- 1.3 ------------------------------------------
% Numerical continuation
% Initializating variables
C_span = 3.09: -0.005 : 3.04;           % Vector of the Jacobi Constant C for numerical contiunation
xx0_halofam = zeros(7,length(C_span));  % Initial contidions for each halo orbit
tf_vec = zeros(1,length(C_span));       % xz crossing time
colors = parula(length(C_span));        % Color variable for the plot

figure;           % Start the figure to plot the finded halo orbits
hold on
grid on
for i = 1:length(C_span)
    C_ref = C_span(i);    % Set the Jacobi constant selected for the i-orbit.
    iter = 0;             % Reset the counter for the number of iteration 
    % Reset the initial tolerances 
    err_y = 1;           
    err_vxf = 1;        
    err_vzf = 1;
    err_C = 1;

    % While cycle to correct the initial conditions
    while (abs(err_vxf) > tol || abs(err_vzf) > tol || abs(err_C) > tol || abs(err_y) > tol) && iter < Nmax
    
        % Evaluate the Jacobi constant and its gradient wrt new initial conditions
        [C_iter,dC_dx0] = JacobiConstant(xx0_new,mu);
    
        % Perform propagation up to event time
        [xf,PHI,te]  = propagate3D(0,xx0_new, tf , mu);
    
        % Evaluate the flow derivative wrt the final time (EOM (dphi/dt)) 
        [dxdt] = xyzCR3BP_STM(0,[xf;PHI(:)], mu);
    
        % Create the transition matrix: MAT [4x4]
        MAT(1:3,1:3) = [PHI(2,1), PHI(2,3) PHI(2,5);           % Add the STM components
                        PHI(4,1), PHI(4,3) PHI(4,5);
                        PHI(6,1), PHI(6,3) PHI(6,5)];
    
        MAT(4,:) =     [dC_dx0(1),dC_dx0(3) dC_dx0(5), 0];     % Add the Jacobin Constant gradient components
    
        MAT(:,4) =     [dxdt(2),dxdt(4),dxdt(6),0];            % Add the time derivative components
    
        % Compute the deviation in the final state
          err_y =   xf(2);
          err_vxf = xf(4);
          err_vzf = xf(6);
          err_C =  C_iter -  C_ref;
          errors = [err_y, err_vxf, err_vzf, err_C]';  % error vector
    
        % Compute the correcction 
          correction = correction - MAT \ errors;
    
        % Update the correction on the initial guess vector
          xx0_new(1) = correction(1);
          xx0_new(3) = correction(2);
          xx0_new(5) = correction(3);
                  tf = correction(4);
    
        % Update iteration counter
        iter = iter+1;
    end

% Get all first guess for the halo orbits varing C
xx0_halofam(:,i) = [xx0_new;tf];
% Propagate forward and backward each orbit and save the result vectors xx
[xf,PHIf,tf, xx_F, tt]  = propagate3D(0,xx0_new,te,mu,false);
[~,~,~, xx_B, ~]  = propagate3D(0,xx0_new,-te,mu,false);

% Plot the propagated orbit
plot3(xx_F(:,1),xx_F(:,2),xx_F(:,3),'LineWidth',1.5,'Color',colors(i,:))
plot3(xx_B(:,1),xx_B(:,2),xx_B(:,3),'LineWidth',1.5,'Color',colors(i,:))
end

% Moon
plot3(P2(1), P2(2), 0, 'Marker', 'o', 'LineWidth', 1, 'Color', 'none', ...
    'MarkerSize', 12, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.7, 0.7, 0.7], 'HandleVisibility', 'off');
text(P2(1), P2(2)-offset/4, 0, 'Moon', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
    'FontWeight', 'bold', 'FontSize', 13);
% L2
plot3(-mu+L2(1), L2(2), 0, 'Marker', 'o', 'LineWidth', 1, 'Color', 'none', ...
    'MarkerSize', 7, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
text(L2(1), L2(2), 0, 'L2', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'FontWeight', 'bold', 'FontSize', 13);
% plot options
axis square
xlabel('x [DU]')
ylabel('y [DU]')
zlabel('z [DU]')
view(20.4915,14.2576)
% Colorbar settings
cbar = colorbar;  
clim([min(C_span), max(C_span)]);  % Limitis
colormap(flipud(colormap)); 
cbar.Label.String = 'Jacobi Constant';  % Title colorbar
title('Halo Orbits varing Jacobi constant $\textit{C}$','FontSize',14)

%% FUNCTIONS

%--------------------------------------------------------------------------
function [dxdt] = xyzCR3BP_STM(~,xx, mu)
% -------------------------------------------------------------------
% Compute the derivative of the state and of th State Transition Matrix xx
%
% INPUT
% xx     [42x1]  State of the system and of th State Transition Matrix
% mu     [1]     Gravitational parameter of the actractor
%
% OUTPUT
% dxdt   [42x1] Derivative of state of the system  and of th State Transition Matrix
%------------------------------------------------------------------------
    
    % Extract variables
    x  = xx(1);
    y  = xx(2);
    z  = xx(3);
    vx = xx(4);
    vy = xx(5);
    % Put PHI in matrix form
    Phi = reshape(xx(7:end),6,6);
    % Compute distances from bodies 1 and 2
    r1 = sqrt((x + mu)^2 + y^2 + z^2);
    r2 = sqrt((x + mu - 1)^2 + y^2 + z^2);
    % Compute derivative of the potential
    dUdx = x - (1-mu)/r1^3*(mu+x) + mu/r2^3*(1-mu-x);
    dUdy = y - (1-mu)/r1^3*y - mu/r2^3*y;
    dUdz =   - (1-mu)/r1^3*z - mu/r2^3*z;
    % Assemble the matrix A(t)=dfdx 6x6 matrix
    dfdx = [                                                                                                                                                                                                                                        0,                                                                                                                                                                                     0,                                                                                                                                                                                 0,  1, 0, 0
                                                                                                                                                                                                                                                    0,                                                                                                                                                                                     0,                                                                                                                                                                                 0,  0, 1, 0
                                                                                                                                                                                                                                                    0,                                                                                                                                                                                     0,                                                                                                                                                                                 0,  0, 0, 1
            (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) + (3*mu*(2*mu + 2*x - 2)*(mu + x - 1))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*(2*mu + 2*x)*(mu + x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2)) + 1,                                                                     (3*mu*y*(mu + x - 1))/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*(mu + x)*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2),                                                                 (3*mu*z*(mu + x - 1))/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*z*(mu + x)*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2),  0, 2, 0
                                                                                                                    (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2)), (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) - (3*y^2*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) + (3*mu*y^2)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) + 1,                                                                                   (3*mu*y*z)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2), -2, 0, 0
                                                                                                                    (3*mu*z*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*z*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2)),                                                                                       (3*mu*y*z)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2), (mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) - (3*z^2*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) + (3*mu*z^2)/((mu + x - 1)^2 + y^2 + z^2)^(5/2),  0, 0, 0];

    

    % Compute the derivative of the STM
    Phidot = dfdx*Phi;
    % Assemble right-hand side
    dxdt = zeros(42,1);

    dxdt(1:3) = xx(4:6);
    dxdt(4)   = dUdx + 2*vy;
    dxdt(5)   = dUdy - 2*vx;
    dxdt(6)   = dUdz;
    dxdt(7:end) = Phidot(:); 
end

%--------------------------------------------------------------------------
function [xf,PHIf,tf, xx, tt]  = propagate3D(t0,x0,tf,mu,varargin)
%------------------------------------------------------------------------
% Propagate an orbit from an initial time t0 to a fineal time tf, given the
% initial state conditions x0 and the mu constant
%
% INPUT
% t0        [1] Initial time
% x0        [6x1] Initial state
% tf        [1] Final time
% mu        [1]  Gravitational parameter of the actractor
% varargin  [1] Event flag   

% OUTPUT
% xf        [6x1] Final state
% PHIf      [6x6] State transition matrix from t0 to tf
% tf        [1] Final time
% xx        [nx6] Propagated state at each time step
% tt        [nx1] Propagation time vector
%--------------------------------------------------------------------------
   % conditons active or not on xz-plane crossing
    if nargin>4
        evtFlag=varargin{1};
    else
        evtFlag=true;
    end
    
    % Time of Flight
    tof = tf - t0;
    % Initialize State Transition Matrix at t0
    Phi0 = eye(6);
    % Append to initial conditions the conditions for the STM
    x0Phi0 = [x0; Phi0(:)];         
    % ODE options
    options_STM = odeset('reltol', 1e-12, 'abstol', 1e-12,'Events',@(x,y) xz_plane_crossing(x,y,evtFlag));
    % ODE solver
    [tt, xx] = ode78(@(t,x) xyzCR3BP_STM(t,x,mu), [0 tof], x0Phi0, options_STM);   % Call xyzCR3BP_STM  

    % Extract state vector and State Transition Matrix
    xf = xx(end,1:6)';
    PHIf = reshape(xx(end,7:end),6,6);     % from vector to STM
    tf = tt(end);

end

%--------------------------------------------------------------------------
function [value, isterminal, direction] = xz_plane_crossing(~,y,isTerminal)
%------------------------------------------------------------------------
% Event function to stop the integration when the xz plane is crossed
%
% INPUT
% y             [6x1] State at a given time step
% isTerminal    [1] Stopping condition
%
% OUTPUT
% value         [1] Second component of the state
% isterminal    [1] Stopping condition
% direction     [1] Direction when stopping criteria is met
%--------------------------------------------------------------------------

% Stop the integration when the xz-plane is reached
    value = y(2);              % Condition on xz plane crossing
    isterminal = isTerminal;   % Terminator
    direction = 0;             % Direction = 0 to consider both crossing sides
end

%--------------------------------------------------------------------------
function [C, dC_dx0] = JacobiConstant(xx,mu)
%--------------------------------------------------------------------------
% From a given state,  evaluate the jacobi constant: 2*OM - v^2
% compute the derivative wrt the given intital state
%
% INPUT
% xx        [6x1] Initial condition
% mu        [1] Gravitational parameter of the main actractor

% OUTPUT
% C         [1x1]  Jacobi constant
% dC_dX0    [1x6] Jacobi constant gradient
%---------------------------------------------------

% Assign value to auxiliary variables
x = xx(1);
y = xx(2);
z = xx(3);
vx = xx(4);
vy = xx(5);
vz = xx(6);

% Potentialc OM evaluation
OM = 1/2*(x^2 + y^2) + (1-mu)/sqrt( (x+mu)^2 +y^2 +z^2) + mu/sqrt( (x+mu-1)^2 +y^2 +z^2) + 1/2*mu*(1-mu);
% Compute Jacobi constant
C = 2*OM - (vx^2 + vy^2 + vz^2);
% Compute the Gradient dJ/dx0
 dC_dx0 = [2*x + ((2*mu + 2*x)*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2) - (mu*(2*mu + 2*x - 2))/((mu + x - 1)^2 + y^2 + z^2)^(3/2), 2*y - (2*mu*y)/((mu + x - 1)^2 + y^2 + z^2)^(3/2) + (2*y*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2), (2*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2) - (2*mu*z)/((mu + x - 1)^2 + y^2 + z^2)^(3/2), -2*vx, -2*vy, -2*vz];
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