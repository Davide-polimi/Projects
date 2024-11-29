function [engineSize, engineSizeCEA] = engineSizing(P_c,T,OF,eps_c,eps_exit,L_star)

% INPUTS in si units:
%       P_c [Pa]
%       T [N]
%       L_star [m]

% OUTPUTS:
    % engineSizing struct with: diameters, length, RAO angles, mass flow rates,
    % specific Impulse
    % VANNO AGGIUNTE LE NOZZLE LOSSES !!!

% Oxidizer & Fuel Data:
temp_RP1 = 298.14;                      % RP-1 Storage Temperature [K]
ent_RP1 = -5430;                        % RP-1 Storage Enthalpy [cal/mol]   
temp_LOX = 90.14;                       % LOX Storage Temperature [K]
ent_LOX = -3032;                        % LOX Storage Enthalpy [cal/mol]

% Constraints/Assumptions
alpha = deg2rad(15);                    % Aperture angle reference conical nozzle [rad]
beta = deg2rad(45);                     % Closure angle reference conical nozzle [rad]

%% Nominal Sizing at t = 0s with nominal pressure
% Solving with CEAM:

engineSizeCEA = CEA('problem','rkt','frozen','o/f',OF,'p,bar',P_c*10^-5,...
    'sup,ae/at',eps_exit,'reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',...
    100,'t,K',temp_RP1,'h,cal/mol',ent_RP1,'oxid','O2(L)','O',2,'wt%',100,'t,K',temp_LOX,...
    'h,cal/mol',ent_LOX,'output', 'short', 'tran','mks','end');

% CEA Outputs:
M_c = 0.066;                                        % Mach number in CC [-] --> from CEA program
engineSizeCEA.output.froz.mach(1) = M_c;
cT = engineSizeCEA.output.froz.cf(3);               % Thrust coefficient [-]
P_e = engineSizeCEA.output.froz.pressure(3)*10^5;   % Exit pressure [bar]
M_e = engineSizeCEA.output.froz.mach(3);
v_e = M_e*engineSizeCEA.output.froz.sonvel(3);
v_t = engineSizeCEA.output.froz.sonvel(2);
rho_t = engineSizeCEA.output.froz.density(2);
mi_t = engineSizeCEA.output.froz.viscosity(2);
k_t = -engineSizeCEA.output.froz.gamma(2)*engineSizeCEA.output.froz.dlvpt(2); 

% Geometrical Sizing:
A_t = T/(cT*P_c);                        % Throat Area [m^2]
D_t = sqrt(4*A_t/pi);                    % Throat Diameter [m]

A_c = eps_c*A_t;                         % Chamber Area [m^2]
D_c = sqrt(4*A_c/pi);                    % Chamber Diameter [m]

A_e = eps_exit*A_t;                      % Exit Area [m^2]
D_e = sqrt(4*A_e/pi);                    % Exit Diameter [m]

V_c = L_star*A_t;                        % Combustion Chamber volume [m^3]

L_c = V_c/A_c;                           % CC Length [m]
L_conv = 0.5*(D_c-D_t)/tan(beta);        % Convergent Length [m]
L_div = 0.7*0.5*(D_e-D_t)/tan(alpha);    % Divergent Length (70% of reference CD nozzle) [m]

% From RAO curves:
th_in = 38;                              % Initial RAO Bell Angle [deg]
th_ex = 9;                               % Exit RAO Bell Angle

lambda = 0.5*(1+cos((alpha+deg2rad(th_ex))/2));

R_curv = 0.382*D_t/2;                       % Radius of curvature of contour nozzle at the throat [m]

Re_t = (R_curv/(D_c/2))^(-1/2)*rho_t*v_t*D_t/mi_t^(-6);
engineSizeCEA.Cd_nozzle = 1-((k_t+1)/2)^(3/4)*(3.266-2.128/(k_t+1))*Re_t^(-1/2)+...
            0.9428*(k_t-1)*(k_t+2)/(k_t+1)^(1/2)*Re_t^(-1);

% Initial Mass Flow Rates:               
mdot = (T-P_e*A_e)/v_e;                  % Propellant Mass Flow Rate [kg/s]
mdot_LOX = mdot*OF/(1+OF);               % Oxidizer Mass Flow Rate [kg/s]
mdot_RP1 = mdot-mdot_LOX;                % Fuel Mass Flow Rate [kg/s]

% Specific Impulse:
Isp = engineSizeCEA.output.froz.isp(3);
Isp_vac = engineSizeCEA.output.froz.isp_vac(3);

% Sizing Storage in 'engineSizing' struct:
engineSize.diameters.Dchamber = D_c;
engineSize.diameters.Dthroat = D_t;
engineSize.diameters.Dexit = D_e;
engineSize.Areas.Achamber = A_c;
engineSize.Areas.Athroat = A_t;
engineSize.Areas.Aexit = A_e;
engineSize.lengths.Lchamber = L_c;
engineSize.lengths.Lconv = L_conv;
engineSize.lengths.Ldiv = L_div;
engineSize.lengths.curvrad = R_curv;
engineSize.RAO.thIn = th_in;
engineSize.RAO.thOut = th_ex;
engineSize.mass.mdot = mdot;
engineSize.mass.mdotLOX = mdot_LOX;
engineSize.mass.mdotRP1 = mdot_RP1;
engineSize.specificImpulse.Isp = Isp;
engineSize.specificImpulse.Isp_vacuum = Isp_vac;
engineSize.twoDlosses.lambda = lambda;



