% Blocco propulsione del loop
clearvars
close all
clc
add_settings

% INPUT
% da LOOP:
% First Stage
variables.stage1.hLaunch = 12000;            % [m]
variables.stage1.RP1LOX.OF = 2.7;            % [-]
variables.stage1.LH2LOX.OF = 6;              % [-]
variables.stage1.METALOX.OF = 3.6;           % [-]
variables.stage1.TWi = 1.3;                  % [-]
variables.stage1.Mp = 17000;                 % [kg]

% Second Stage
variables.stage2.hLaunch = 80e3;            % [m]
variables.stage2.RP1LOX.OF = 2.7;            % [-]
variables.stage2.LH2LOX.OF = 6;              % [-]
variables.stage2.METALOX.OF = 3.6;           % [-]
variables.stage2.TWi = 0.95;                 % [-]
variables.stage2.Mp = 2000;                  % [kg]

% da BASELINE:
%   prop couple
%   Pcc
%   eps_UpStage
%   MStack

% Stage 1
% Pressure in combustion chamber
baseline.stage1.Pcc_bar = 70;              % [bar]
% GLOM al passo k-1
baseline.stage1.Mstack = 21500;            % [kg]

% Stage 2
baseline.stage2.Pcc_bar = 50;              % [bar]             
baseline.stage2.Mstack = 3000;             % [kg]
baseline.stage2.epsUS = 150;               

% Target Orbit
baseline.hOrbit = 400e3;                    % [m]


% DA DEFINIRE COME PASSARLO DENTRO AL CODICE

% Nella Baseline devo inserire tutte le variabili che servono per il CEA


% OUTPUT
%   Isph
%   Th
%   mdot_p
%   tburn
%   Ptank
%   Geom nozzle
[propulsion] = propSystem(variables,baseline,{'1'},{'METALOX'});

%% FUNCTION

function [propulsion] = propSystem(variables,baseline,stages,propCouple)
% PropSystems computes the key parameter of the propulsion
% 
% INPUT:
% variables         : Structure containing all the variables of the loop
% baseline          : Structure containing all the parameters coming from
%                     baseline
% stages            : Cell array of strings specifying the stages to compute
%                     (e.g., {'1', '2', '3'})
% propCouple        : Cell array of strings specifying the propellant type for each stage
%                     (e.g., {'RP1-LOX', 'LH2-LOX', 'METALOX'})% 
% OUTPUT:
% propulsion        : Structure containing all the outputs needed from Prop Block
% 

% Initialize Propulsion as an empty structure
propulsion = struct();

% Check length consistency
if length(stages) ~= length(propCouple)
    error('The number of stages and propellant types must be the same.');
end

% Loop through each stage
for i = 1:length(stages)
    nStage = stages{i};
    couple = propCouple{i};

switch nStage
        case {'first','1'}
            switch couple
                case {'RP1-LOX','LH2-LOX','METALOX'}
                    propulsion.stage1 = calculateFirstStageLiquid(variables,baseline,couple);
                case 'solid'
                    propulsion.stage1 = calculateFirstStageSolid(variables,baseline);
                otherwise
                    error('Wrong propCouple for First Stage. Use "RP1-LOX", "LH2-LOX", "METALOX" or "solid".');
            end

        case {'second','2'}
            switch couple
                case {'RP1-LOX','LH2-LOX','METALOX'}
                    propulsion.stage2 = calculateSecondStageLiquid(variables,baseline,couple);
                case 'solid'
                    propulsion.stage2 = calculateSecondStageSolid(variables,baseline);    
                otherwise
                    error('Wrong propCouple for Second Stage. Use "RP1-LOX", "LH2-LOX", "METALOX" or "solid".');
            end
            
        case {'third','3'}
            switch couple
                case 'RP1-LOX'
                    disp('Third Stage with RP1-LOX propellant.');
                    % Codice per il primo stadio, motore liquido
                    propulsion.stage3 = calculateThirdStageLiquid(variables, baseline);
                case 'LH2-LOX'  
                     disp('Third Stage with LH2-LOX propellant.');
                    % Codice per il primo stadio, motore liquido
                    propulsion.stage3 = calculateThirdStageLiquid(variables, baseline);
                case 'METALOX'
                     disp('Third Stage with METALOX propellant.');
                    % Codice per il primo stadio, motore liquido
                    propulsion.stage3 = calculateThirdStageLiquid(variables, baseline);

                otherwise
                    error('Wrong propCouple for Third Stage. Use "RP1-LOX", "LH2-LOX", "METALOX".');
            end
            
otherwise
            error('Wrong nStage. Use "first" or "1", "second" or "2", "third" or "3".');
end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stage1 = calculateFirstStageLiquid(variables,baseline,couple)

% Extract the variables
hLaunch = variables.stage1.hLaunch;                 % [km]
TW_i = variables.stage1.TWi;                        % [-]
Mp = variables.stage1.Mp;                           % [kg]
hEnd = variables.stage2.hLaunch;                    % detach First stage [km]

% Pressure at the exit of the nozzle
[~,~,Pe,~,~] = atmosisa(hLaunch,'extended','on');   % [Pa]
Pe_bar = Pe * 1e-5;                                 % [bar]

% Extract parameters from baseline
Pcc_bar = baseline.stage1.Pcc_bar;          % Pressure in combustion chamber [bar]
Pcc = Pcc_bar * 10^5;                       % Pressure in combustion chamber [Pascal]
MStack = baseline.stage1.Mstack;            % Mass of the stack [kg]

% Efficiency
etaCstar = 0.98;
etaCt = 0.98;
etaIsp = etaCstar * etaCt;

% Constants
g0 = 9.8065;        % Gravity
Ru = 8314.46;       % Universal Gas Constant

switch couple
       case 'RP1-LOX'
                    disp('First Stage with RP1-LOX propellant.');
                    % Extract the correspondant OF Ratio
                    OF = variables.stage1.RP1LOX.OF;           % [-]
                    
                    % Through CEA compute all the thermophysical properties
                    outCEAM = CEA('prob','hp','p(bar)',Pcc_bar,'o/f',OF,...
                                  'reac','fuel','RP-1(L)','C',1.,'H',1.9423,'wt%',100, 'h,cal/mol',-5430.,'t(k)',300.0, ...
                                  'oxid','O2(L)',...
                                  'output','short','transport','end');
                    
                    % The characteristic length Lstar depends on the
                    % propellant couple chosen
                    Lstar = 1.143;              % [m]

       case 'LH2-LOX'  
                     disp('First Stage with LH2-LOX propellant.');
                    % Extract the correspondant OF Ratio
                    OF = variables.stage1.LH2LOX.OF;           % [-]
                    
                    % Through CEA compute all the thermophysical properties
                    outCEAM = CEA('prob','hp','p(bar)',Pcc_bar,'o/f',OF,...
                                  'reac','fuel','H2(L)', ...
                                  'oxid','O2(L)',...
                                  'output','short','transport','end');

                    % The characteristic length Lstar depends on the
                    % propellant couple chosen
                    Lstar = 0.9;              % [m]
       case 'METALOX'
                    disp('First Stage with METALOX propellant.');
                    % Extract the correspondant OF Ratio
                    OF = variables.stage1.METALOX.OF;           % [-]
                    
                    % Through CEA compute all the thermophysical properties
                    outCEAM = CEA('prob','hp','p(bar)',Pcc_bar,'o/f',OF,...
                                  'reac','fuel','CH4(L)', 'C',1.,'H',4, 'wt%', 100,...
                                  'oxid','O2(L)',...
                                  'output','short','transport','end');


                    % The characteristic length Lstar depends on the
                    % propellant couple chosen
                    Lstar = 1;              % [m]
end

% Extract the parameters from CEA output
Tcc = outCEAM.output.temperature;
Mmol = outCEAM.output.mw;
Cp = outCEAM.output.cp_tran.froz*1000;
    
R = Ru/Mmol;                     % Costante particolare del gas
gamma = Cp/(Cp-R);               % Specific Heat Ratio

% Nozzle Exit Velocity
ve = sqrt( 2*gamma/(gamma-1)* R * Tcc * (1 - (Pe_bar/Pcc_bar)^((gamma-1)/gamma)) );

% Nozzle Exit Temperature considering adiabatic expansion model 
Te = Tcc * (Pe_bar/Pcc_bar) ^ (gamma-1)/gamma;

% Nozzle Exit Mach
ae = sqrt(gamma*R*Te);           % [m/s]
Me = ve/ae;                      % [-]

% Thrust in optimal expansion condition
T_opt = TW_i * MStack * g0;      % [N]

% Since we are in optimal expansion cond. at the detachment from aircraft,
% I can recover the mass flow rate (assumed constant)
mdot_p = T_opt/ve;               % [kg/s]

% Calcolo Cstar ideale con i dati del CEA, e conoscendo etaCstar mi trovo
% il valore vero
Cstar = sqrt((R*Tcc)/(gamma*((2/(gamma+1))^((gamma+1)/(gamma-1)))));
Cstar = Cstar * etaCstar;

% Specific Impulse in optimal expansion condition
Isp = T_opt / (mdot_p*g0);
Isp_opt = Isp * etaIsp;

% Burning Time
tburn = Mp/mdot_p;

% Conoscendo l'h del launcher, posso calcolare la variazione di T e Isp
% durante la missione
hRange = linspace(hLaunch,hEnd,10000);
Pamb = zeros(length(hRange),1);

for i = 1:length(hRange)
    [~,~,Pamb(i),~,~] = atmosisa(hRange(i),"extended","on",'action','None');

end

% Nozzle
eps = ( ((gamma+1)/2)^(1/(gamma+1)) * (Pe_bar/Pcc_bar)^(1/gamma) * sqrt( (gamma+1)/(gamma-1)*(1 - (Pe_bar/Pcc_bar)^((gamma-1)/gamma))) )^-1;

At = Cstar * mdot_p / Pcc;
Ae = eps * At;

% thrust e Isp in funzione dell'altitudine
T = T_opt + (Pe - Pamb)*Ae;
Isp = Isp_opt + (Pe - Pamb) * Ae/(mdot_p * g0);

% Plot the evolution of Thrust and Isp over the altitude
figure()
subplot(1,2,1)
plot(hRange/10^3,T.*ones(length(hRange),1))
title('Thrust Profile','FontSize',20)
xlabel('Altitude [km]')
ylabel('Thrust [N]','Rotation',0)

subplot(1,2,2)
plot(hRange/10^3,Isp.*ones(length(hRange),1),'Color','g')
title('Isp Profile','FontSize',20)
xlabel('Altitude [km]')
ylabel('Isp [s]','Rotation',0)

sgtitle('I Stage');

% GEOMETRY
alphaDiv = deg2rad(15);     % Impose alpha = 15° [rad]
betaDiv = deg2rad(45);      % Impose beta = 45° [rad]
Mcc = 0.3;                  % Mach in CC [-]

% COMBUSTION CHAMBER
Acc = At/Mcc * ( (2/(gamma+1)) * (1-((gamma-1)/2 * Mcc^2)));    % Area [m^2]
% Guess of the Volume of the combustion chamber using the Lstar
Vcc = At * Lstar;               % Volume [m^3]
Lcc = Vcc/Acc;                  % Length [m]

% NOZZLE 
Dcc = sqrt(Acc/pi * 4);     % Diameter CC [m]
Dt = sqrt(At/pi * 4);       % Throat Diameter [m]
De = sqrt(Ae/pi * 4);       % Nozzle Exit Diameter [m]

% Assume the convergent part of the bell nozzle is the same as the
% reference conical one
Lconv = 0.5 * (Dcc - Dt)/tan(betaDiv);     % Length of the convergent part of the nozzle [m]
Ldiv_conical = 0.5 * (De - Dt)/tan(alphaDiv);
Ltot = Lcc + Lconv + Ldiv_conical;

% OUTPUT
stage1.Me = Me;                 % Mach at the exit of the nozzle
stage1.T_opt = T_opt;           % Thrust in optimal expansion condition
stage1.T = T;                   % Vector containing thrust at each altitude
stage1.mdot_p = mdot_p;         % Mass flow rate
stage1.Isp_opt = Isp_opt;       % Isp in optimal expansion condition
stage1.Isp = Isp;               % Vector containing Isp at each altitude
stage1.hRange = hRange;         % Vector containing the altitude 
stage1.tburn = tburn;           % Burning Time
stage1.geometry.Dt = Dt;
stage1.geometry.De = De;
stage1.geometry.Dcc = Dcc;
stage1.geometry.Lcc = Lcc;
stage1.geometry.Lconv = Lconv;
stage1.geometry.Ldiv = Ldiv_conical;
stage1.geometry.Ltot = Ltot;
stage1.eps = eps;
end

% ---------------------------------------------------------------------------------------------------- %
% function [stage1] = calculateFirstStageSolid(variables,baseline)



% ---------------------------------------------------------------------------------------------------- %
function [stage2] = calculateSecondStageLiquid(variables,baseline,couple)

% Extract the variables
TW_i = variables.stage2.TWi;               % [-]
Mp = variables.stage2.Mp;                  % [kg]
hLaunch = variables.stage2.hLaunch;        % Altitude of stage 1 detachment [m]

% Pressure at the exit of the nozzle
[~,~,Pe,~,~] = atmosisa(hLaunch);   % [Pa]
Pe_bar = Pe * 1e-5;                 % [bar]

% Extract parameters from baseline
Pcc_bar = baseline.stage2.Pcc_bar;         % Pressure in combustion chamber [bar]
Pcc = Pcc_bar * 10^5;                      % Pressure in combustion chamber [Pascal]
MStack = baseline.stage2.Mstack;           % Mass of the stack [kg]

hOrbit = baseline.hOrbit;
% Efficiency
etaCstar = 0.98;
etaCt = 0.98;
etaIsp = etaCstar * etaCt;

% Constants
g0 = 9.8065;        % Gravity
Ru = 8314.46;       % Universal Gas Constant

switch couple
       case 'RP1-LOX'
                    disp('Second Stage with RP1-LOX propellant.');
                    % Extract the correspondant OF Ratio
                    OF = variables.stage2.RP1LOX.OF;           % [-]
                    
                    % Through CEA compute all the thermophysical properties
                    outCEAM = CEA('prob','hp','p(bar)',Pcc_bar,'o/f',OF,...
                                  'reac','fuel','RP-1(L)','C',1.,'H',1.9423,'wt%',100, 'h,cal/mol',-5430.,'t(k)',300.0, ...
                                  'oxid','O2(L)',...
                                  'output','short','transport','end');
                    
                    % The characteristic length Lstar depends on the
                    % propellant couple chosen
                    Lstar = 1.143;              % [m]

       case 'LH2-LOX'  
                     disp('Second Stage with LH2-LOX propellant.');
                    % Extract the correspondant OF Ratio
                    OF = variables.stage2.LH2LOX.OF;           % [-]
                    
                    % Through CEA compute all the thermophysical properties
                    outCEAM = CEA('prob','hp','p(bar)',Pcc_bar,'o/f',OF,...
                                  'reac','fuel','H2(L)', ...
                                  'oxid','O2(L)',...
                                  'output','short','transport','end');

                    % The characteristic length Lstar depends on the
                    % propellant couple chosen
                    Lstar = 0.9;              % [m]
       case 'METALOX'
                    disp('Second Stage with METALOX propellant.');
                    % Extract the correspondant OF Ratio
                    OF = variables.stage2.METALOX.OF;           % [-]
                    
                    % Through CEA compute all the thermophysical properties
                    outCEAM = CEA('prob','hp','p(bar)',Pcc_bar,'o/f',OF,...
                                  'reac','fuel','CH4(L)', 'C',1.,'H',4, 'wt%', 100,...
                                  'oxid','O2(L)',...
                                  'output','short','transport','end');


                    % The characteristic length Lstar depends on the
                    % propellant couple chosen
                    Lstar = 1;              % [m]
end

% Extract the parameters from CEA output
Tcc = outCEAM.output.temperature;
Mmol = outCEAM.output.mw;
Cp = outCEAM.output.cp_tran.froz*1000;
    
R = Ru/Mmol;                     % Costante particolare del gas
gamma = Cp/(Cp-R);               % Specific Heat Ratio

% Nozzle Exit Velocity
ve = sqrt( 2*gamma/(gamma-1)* R * Tcc * (1 - (Pe_bar/Pcc_bar)^((gamma-1)/gamma)) );

% Nozzle Exit Temperature considering adiabatic expansion model 
Te = Tcc * (Pe_bar/Pcc_bar) ^ (gamma-1)/gamma;

% Nozzle Exit Mach
ae = sqrt(gamma*R*Te);           % [m/s]
Me = ve/ae;                      % [-]

% Thrust in optimal expansion condition
T_opt = TW_i * MStack * g0;      % [N]

% Since we are in optimal expansion cond. at the detachment from aircraft,
% I can recover the mass flow rate (assumed constant)
mdot_p = T_opt/ve;               % [kg/s]

% Calcolo Cstar ideale con i dati del CEA, e conoscendo etaCstar mi trovo
% il valore vero
Cstar = sqrt((R*Tcc)/(gamma*((2/(gamma+1))^((gamma+1)/(gamma-1)))));
Cstar = Cstar * etaCstar;

% Specific Impulse in optimal expansion condition
Isp = T_opt / (mdot_p*g0);
Isp_opt = Isp * etaIsp;

% Burning Time
tburn = Mp/mdot_p;

% Nozzle
eps = ( ((gamma+1)/2)^(1/(gamma+1)) * (Pe_bar/Pcc_bar)^(1/gamma) * sqrt( (gamma+1)/(gamma-1)*(1 - (Pe_bar/Pcc_bar)^((gamma-1)/gamma))) )^-1;

if eps > baseline.stage2.epsUS
   eps  = baseline.stage2.epsUS;
end

% GEOMETRY
At = Cstar * mdot_p / Pcc;
Ae = eps * At;

% For the second stage, due to the altitude reached, pamb is negligible
T = T_opt + Pe * Ae;
Isp = Isp_opt + Pe * Ae/(mdot_p * g0);

% GEOMETRY
alphaDiv = deg2rad(15);     % Impose alpha = 15° [rad]
betaDiv = deg2rad(45);      % Impose beta = 45° [rad]
Mcc = 0.3;                  % Mach in CC [-]

% COMBUSTION CHAMBER
Acc = At/Mcc * ( (2/(gamma+1)) * (1-((gamma-1)/2 * Mcc^2)));    % Area [m^2]
% Guess of the Volume of the combustion chamber using the Lstar
Vcc = At * Lstar;               % Volume [m^3]
Lcc = Vcc/Acc;                  % Length [m]

% NOZZLE 
Dcc = sqrt(Acc/pi * 4);     % Diameter CC [m]
Dt = sqrt(At/pi * 4);       % Throat Diameter [m]
De = sqrt(Ae/pi * 4);       % Nozzle Exit Diameter [m]

% Assume the convergent part of the bell nozzle is the same as the
% reference conical one
Lconv = 0.5 * (Dcc - Dt)/tan(betaDiv);     % Length of the convergent part of the nozzle [m]
Ldiv_conical = 0.5 * (De - Dt)/tan(alphaDiv);
Ltot = Lcc + Lconv + Ldiv_conical;

% % Plot the evolution of Thrust and Isp over the altitude
% hRange = linspace(hLaunch,hOrbit,10000);
% 
% figure()
% subplot(1,2,1)
% plot(hRange/10^3,T)
% title('Thrust Profile','FontSize',20)
% xlabel('Altitude [km]')
% ylabel('Thrust [N]','Rotation',0)
% 
% subplot(1,2,2)
% plot(hRange/10^3,Isp,'Color','g')
% title('Isp Profile','FontSize',20)
% xlabel('Altitude [km]')
% ylabel('Isp [s]','Rotation',0)
% 
% sgtitle('II Stage');

% OUTPUT
stage2.Me = Me;                 % Mach at the exit of the nozzle
stage2.T_opt = T_opt;           % Thrust in optimal expansion condition
stage2.T = T;                   % Vector containing thrust at each altitude
stage2.mdot_p = mdot_p;         % Mass flow rate
stage2.Isp_opt = Isp_opt;       % Isp in optimal expansion condition
stage2.Isp = Isp;               % Vector containing Isp at each altitude
stage2.tburn = tburn;           % Burning Time
stage2.geometry.Dt = Dt;
stage2.geometry.De = De;
stage2.geometry.Dcc = Dcc;
stage2.geometry.Lcc = Lcc;
stage2.geometry.Lconv = Lconv;
stage2.geometry.Ldiv = Ldiv_conical;
stage2.geometry.Ltot = Ltot;
stage2.eps = eps;
end

% ---------------------------------------------------------------------------------------------------- %
% function [stage2] = calculateSecondStageSolid(variables,baseline)



% ---------------------------------------------------------------------------------------------------- %
% function [stage3] = calculateThirdStage(variables,baseline,couple)
% 
% 
% 
% 
% end
