function [L,D] = Aerodynamics_computation(rocket_design,v,alpha,h,engine_mode,phi,wing_design)
    % Inputs:
    % 1.Rocket_design
    % - .a,.b: semimajor and semiminor axis of cross section of the body [m]
    % - .x: coordinate at which cross section (a,b) are provided
    % - .x(1) is nose tip, .x(end) is lenght of the body
    % - .dnose: nose diameter
    % - .lnose: nose length
    % - .nosetype: 'char' - ogive,cone...
    % 2. v: velocity
    % 3. alpha:angle of attack %[rad]
    % 4. altitude
    % 5. engine_mode: 1 engine active; 0 coast phase (for base drag evaluation)
    % 6. phi: rotation of the body wrt its axis
    % 7. wing_design: D,shapes of aerodynamic surfaces

    [~,a_sound,~,rho,~] = atmosisa(h,extended="on");
    M = v/a_sound;
    [Ca,~,~,~] = Ca_computation(rocket_design.a,rocket_design.b,rocket_design.x,rocket_design.ln,rocket_design.d_nose,rocket_design.nose_type,M,alpha,h,engine_mode);
    Cn = Cn_computation(rocket_design.a,rocket_design.b,rocket_design.x,rocket_design.ln,rocket_design.d_nose,rocket_design.nose_type,M,alpha,h,phi);
    [Cl,Cd] = aero_conversion(Cn,Ca,alpha);
    d = max(rocket_design.a);
    S = pi*d^2/4;
    L = 1/2*rho*v^2*S*Cl;
    D = 1/2*rho*v^2*S*Cd;
end






function [Ca,Ca_w,Ca_f,Ca_b] = Ca_computation(a,b,x,ln,d_nose,nose_type,M,alpha,h,engine_mode)
% This function computes preliminary aerodynamic axial coefficients
% according to the concept of component build up



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
end

function Cn = Cn_computation(a,b,x,ln,d_nose,nose_type,M,alpha,h,phi,delta,A_wing,Asp_ratio_w)
d = max(a); % reference diameter assigned as the maximum semimajor axis
l = x(end);
A_ref = pi*d^2/4;

if nargin < 10
    phi = 0; % assume zero inclination of elliptical shape if not given
end

if nargin < 11
    Cn_w = 0;
end

% Cn body
Cn_ratio_SB = @(a,b) a/b*cos(phi)^2 + b/a*sin(phi)^2;
Cn_ratio_N = @(a,b) (3/2) * sqrt(a/b) * ((-b^2 / a^2) / (1 - (b^2 / a^2)^3) * log((a / b) * (1 + sqrt(1 - (b^2 / a^2))) + ...
        1 / (1 - (b^2 / a^2))));

n = max(size(a));
r = sqrt(a.*b);

dA = pi*(a(2)*b(2) - a(1)*b(1));
int_SB = (Cn_ratio_SB(a(2),b(2)) + 1) * dA/2;
int_N = Cn_ratio_SB(a(2),b(2))*r(2)*(x(2)-x(1))/2;



for i = 3:n
    dA = pi*(a(i)*b(i) - a(i-1)*b(i-1));
    int_SB = int_SB + (Cn_ratio_SB(a(i),b(i)) + Cn_ratio_SB(a(i-1),b(i-1))) * dA/2;
    if a~=b
        int_N = int_N + (Cn_ratio_N(a(i),b(i))*r(i) + Cn_ratio_N(a(i-1),b(i-1))*r(i-1))*(x(i)-x(i-1))/2;  
    end
end

if M <= 1
    eta = 0.05*(l/d)+0.52;
elseif M > 1
    eta = 1;
end

M_n = M*sin(alpha);
[~,a,~,~,nu] = atmosisa(h,extended="on");
v = M*a;
Re = v*d/nu;
Re_n = Re*sin(alpha);

%%%%%%%%%%%%%%
% ADD correlation for critical reynolds !!!!!!!!!!!!!!!!
%%%%%%%%%%%%%%


% Graphical interpolation of crossflow drag coefficient as a function of
% crossflow reynolds

Cdn_Mn = [0, 1.217;
        0.1, 1.217;
        0.2, 1.217;
        0.3, 1.223;
        0.4, 1.282;
        0.5, 1.363;
        0.6, 1.487;
        0.7, 1.537;
        0.8, 1.446;
        0.9, 1.747;
        1.0, 2.051;
        1.1, 1.819;
        1.2, 1.687;
        1.3, 1.585;
        1.4, 1.569;
        1.5, 1.545;
        1.6, 1.511;
        1.7, 1.493;
        1.8, 1.463;
        1.9, 1.439;
        2.0, 1.432;
        2.1, 1.423;
        2.2, 1.413;
        2.3, 1.399;
        2.4, 1.388;
        2.5, 1.385;
        2.6, 1.382;
        2.7, 1.375;
        2.8, 1.368;
        2.9, 1.368;
        3.0, 1.367;
        3.1, 1.362;
        3.2, 1.357;
        3.3, 1.351;
        3.4, 1.346;
        3.5, 1.341;
        3.6, 1.337;
        3.7, 1.337;
        3.8, 1.336;
        3.9, 1.330;
        4.0, 1.326;
        4.1, 1.326;
        4.2, 1.324;
        4.3, 1.317;
        4.4, 1.311;
        4.5, 1.311];

if M_n <= 4.5
    % Interpolate to find the corresponding C_dn value
    C_dn = interp1(Cdn_Mn(:,1), Cdn_Mn(:,2), M_n, 'linear');
elseif M_n > 4.5
    C_pstag = 1.8;
    C_dn = 2/3*C_pstag;
end

Cn_b = (sin(2*alpha)*cos(alpha/2)*int_SB + 2*eta*C_dn*sin(alpha)^2*int_N)/A_ref;

% Cn wing
if nargin == 13
    alpha_att = alpha + delta; % wing effective angle of attack
    if M^2 >= 1+(8/pi/Asp_ratio_w)^2
        Cn_w = (4*sin(alpha_att)*cos(alpha_att)/(M^2-1)^1/2+2*sin(alpha_att)^2)*A_wing/A_ref;
    elseif M^2 < 1+(8/pi/Asp_ratio_w)^2
        Cn_w = ((pi*Asp_ratio_w/2)*sin(alpha_att)*cos(alpha_att)+2*sin(alpha_att)^2)*A_wing/A_ref;
    end
end

Cn = Cn_b + Cn_w;

end



function [Cl,Cd] = aero_conversion(Cn,Ca,alpha)
    % To obtain lift and drag coefficients from normal and axial coefficients

    Cl = Cn*cos(alpha)-Ca*sin(alpha);
    Cd = Cn*sin(alpha)+Ca*cos(alpha);
end
function cdw = wvdrgogv(diameter,length,Mach)

% Returns the supersonic wave drag for a tangent ogive
% by empirical method of Miles
%
% diameter = ogive base diameter
% length = lengthof ogive
% Mach = Mach number

% The values resulting from these calculations are valid in the 1.5 to 3.5
% Mach range and for semi-vertex angles between 10 and 25 degrees

if Mach > 3.5
    Mach = 3.5; % modifying Mach number into the range of model validity
end

sigma = 2*(180/pi)*atan(diameter/2/length);
P = (0.083+0.096/Mach^2)*(sigma/10)^1.69;
lod2 = (length/diameter)^2;
cdw = P*(1-(196*lod2-16)/(14*(Mach+18)*lod2));
end