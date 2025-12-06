% Launch project Aerodynamics
function [Ca,Ca_w,Ca_f,Ca_b] = Ca_computation(a,b,x,ln,d_nose,nose_type,M,alpha,h,engine_mode)
% This function computes preliminary aerodynamic axial coefficients
% according to the concept of component build up

% Inputs:
% - a,b: semimajor and semiminor axis of cross section of the body [m]
% - x: coordinate at which cross section (a,b) are provided
% x(1) is nose tip, x(end) is lenght of the body
% - dnose
% - lnose
% - nosetype: 'char' - ogive,cone...
% - M
% - alpha % [rad]
% - altitude
% - engine_mode: 1 engine active; 0 coast phase (for base drag evaluation)

%%%%%%%%%%%%%%%%%
% POSSIBLE MODIFICATION, USE INTERPOLATION FOR 1<M<1.3
%%%%%%%%%%%%%%%

if nargin < 10
    engine_mode = 1; % assume engine active if not given as an input
end

if h > 84000
    Ca = 0; % After this altitude the model consider not present Lift and Drag
    Ca_w = 0;
    Ca_f = 0;
    Ca_b = 0;
else
    l = x(end);
    
    fn_nose = ln/d_nose; % fineness nose ratio
    
    % Angle of attack checks
    if alpha <= deg2rad(90) && alpha >= 0
        alpha = + alpha;
    elseif alpha <= deg2rad(180) && alpha >= deg2rad(90)
        alpha = deg2rad(180) - alpha;
    end
    
    if alpha <= deg2rad(90)
        theta = atan(1/2/fn_nose); % cone-half angle
    elseif alpha > deg2rad(90)
        theta = pi/2;
    end
    
    % ------ 1. WAVE DRAG -------
    
    switch nose_type
        case 'C' % conical
            if M <= 1
                Ca_w = 0.8*sin(theta)^2;
    
            elseif M > 1
                beta = sqrt(M^2-1);
    
                Ca_w = (4*sin(theta))^2*(2.5+8*beta*sin(theta))/(1+16*beta*sin(theta)); % Linnell-Bailey experimental relation
            end
        
        case 'TO' % tanget ogive
            if M <= 1
                Ca_w = 0.8*sin(theta)^2;
            
            elseif M > 1
                Ca_w = wvdrgogv(d_nose,ln,M); % function that computes wave drag of ogive nose
            end    
    end
    
    % ------ 2. BASE PRESSURE ---------
    gamma = 1.4;
    if M >= 1
        Cp_b = 2/(gamma*M^2)*((2/(gamma+1))^1.4*(1/M)^2.8*(2*gamma*M^2-(gamma-1))/(gamma+1)-1); % Gabeaud base pressure for blunt based bodies in supersonic flow
        Ca_b = -Cp_b;
    elseif M < 1
        Ca_b = 0.12 + 0.13*M^2; % Fleeman relation for subsonic condition
    end

    if engine_mode == 1
        Cp_b = 0;
    end
    
    % ----- 3. SKIN-FRICTION -------
    Cf_i = 1.48e-2; % Re < 10^4 OpenRocket
    
    % Cf_i = 0.0032*(R_s/l)^0.2; % Re < 10^4 < Re < Re_crit (10^5)
    
    if M >= 1
        C_sf = Cf_i/((1+0.144*M^2)^0.65);
    elseif M<1
        C_sf = Cf_i/(1+0.08*M^2);
    end
    
    A_wet_body = (l-ln)*d_nose + (ln*d_nose)/2; % Area wetted by the flow considering a conical nose
    A_ref = pi*(d_nose/2)^2;
    
    Ca_f = C_sf*(1+1/2/(l/d_nose)*A_wet_body)/A_ref;
    
    % ----- 4. SUM OF CONTRIBUTIONS -------
    Ca = Ca_w + Ca_f + Ca_b;
    Ca = Ca*1.1; % to take into account other drag sources (like parasitic drag)
    
    Ca = Ca*cos(alpha)^2;
end







