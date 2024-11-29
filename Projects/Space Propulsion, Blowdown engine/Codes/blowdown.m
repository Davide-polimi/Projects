function [bd] = blowdown(engineSize,engineSizeCEA,t_end,dt,K_inj_ox, ...
                        K_idr_ox,K_dyn_ox,K_inj_fuel,K_idr_fuel,K_dyn_fuel, ...
                        T_in,P_c,P_c_min,eps_exit,temp_RP1,temp_LOX,ent_RP1, ...
                        ent_LOX,rho_fuel,rho_ox,R,gamma_He,tank,cooling, ...
                        coolingSize)

% Inizialization

bd.P_cc = [];              % Combustion chamber pressure
bd.O_F = [];               % Ox/fuel ratio
bd.T = [];                 % Thrust
bd.Ptank_ox_v = [];        % Oxidizer tank pressure
bd.Ptank_fuel_v = [];      % Fuel tank pressure
bd.V_ox = [];              % Oxidizer tank volume
bd.V_fuel = [];            % Fuel tank volume
bd.Vgas_ox  = [];          % Gas volume in oxidizer tank
bd.Vgas_fuel = [];         % Gas volume in fuel tank
bd.mdot_v = [];            % Mass flow rate
bd.mdot_ox_v = [];         % Oxidizer mass flow rate
bd.mdot_fuel_v = [];       % Fuel mass flow rate
bd.T_cc = [];              % Combustion chamber temperature
bd.m_mol_v = [];           % Molar mass
bd.cstar_v = [];           % Characteristic velocity
bd.cT_v = [];              % Thrust coefficient
bd.k_v = [];               % Propellant gamma
bd.cooling.T_end_cooling = [];
bd.cooling.Delta_P_cool = [];
bd.t_burn = [];
bd.P_e = [];

bd.T(1) = T_in;
bd.P_cc(1) = P_c;
bd.O_F(1) = engineSizeCEA.output.oxfl;
bd.mdot_v(1) = engineSize.mass.mdot;
bd.mdot_ox_v(1) = engineSize.mass.mdotLOX;          
bd.mdot_fuel_v(1) = engineSize.mass.mdotRP1;         
bd.Ptank_ox_v(1) = tank.P_ox;         
bd.Ptank_fuel_v(1) = tank.P_fuel;       
bd.V_ox(1) = tank.V_ox_in;
bd.V_fuel(1) = tank.V_fuel_in;
bd.Vgas_ox(1) = tank.Vgas_ox_in;
bd.Vgas_fuel(1) = tank.Vgas_fuel_in;
bd.k_v(1) = -engineSizeCEA.output.froz.gamma(1)*engineSizeCEA.output.froz.dlvpt(1);     
bd.T_cc(1) = engineSizeCEA.output.froz.temperature(1);                    
bd.m_mol_v(1) = engineSizeCEA.output.froz.mw(1)*10^(-3);                
bd.cT_v(1) = engineSizeCEA.output.froz.cf(3);
bd.cstar_v(1) = sqrt(R/bd.m_mol_v(1)*bd.T_cc(1)/(bd.k_v(1)*(2/(bd.k_v(1)+1))^((bd.k_v(1)+1)/(bd.k_v(1)-1))));
bd.cooling.T_end_cooling(1) = coolingSize.temperatures.T_endChamber;
bd.cooling.Delta_P_cool(1) = coolingSize.losses.DP_cool;
bd.t_burn(1) = 0;
bd.P_e(1) = engineSizeCEA.output.froz.pressure(3)*10^5;

t_simulation = 0:dt:t_end;

for i = 1:length(t_simulation)
% CEAM
p_cycle = CEA('problem','rkt','frozen','o/f',bd.O_F(i),'p,bar',bd.P_cc(i)*10^-5,...
'sup,ae/at',eps_exit,'reactants','fuel','RP-1(L)','C',1,'H',1.9423,'wt%',...
100,'t,K',temp_RP1,'h,cal/mol',ent_RP1,'oxid','O2(L)','O',2,'wt%',100,'t,K',temp_LOX,...
'h,cal/mol',ent_LOX,'output', 'short', 'tran','mks','end');

bd.T_cc(i) = p_cycle.output.froz.temperature(1);
bd.k_v(i) = -p_cycle.output.froz.gamma(1)*p_cycle.output.froz.dlvpt(1);
bd.cT_v(i) = p_cycle.output.froz.cf(3);
bd.m_mol_v(i) = p_cycle.output.froz.mw(1)*10^-3;
bd.cstar_v(i) = sqrt(R/bd.m_mol_v(i)*bd.T_cc(i)/(bd.k_v(i)*(2/(bd.k_v(i)+1))^((bd.k_v(i)+1)/(bd.k_v(i)-1))));
bd.P_e(i) = p_cycle.output.froz.pressure(3)*1e5;
% Volume tanks
bd.V_ox(i+1) = bd.V_ox(i)-bd.mdot_ox_v(i)/rho_ox*dt;
bd.V_fuel(i+1) = bd.V_fuel(i)-bd.mdot_fuel_v(i)/rho_fuel*dt;

if bd.P_cc(i) <= P_c_min*1.02
    bd.t_burn(i) = 0;
elseif bd.V_ox(i+1) <= tank.V_ox_in*0.1 || bd.V_fuel(i+1) <= tank.V_fuel_in*0.1
    bd.t_burn(i+1) = 0;
else
    bd.t_burn(i) = i*dt;
end

% Gas volume
bd.Vgas_ox(i+1) = tank.V_ox_tot-bd.V_ox(i+1);
bd.Vgas_fuel(i+1) = tank.V_fuel_tot-bd.V_fuel(i+1);

% Pressure tanks (Adiabatic)
bd.Ptank_ox_v(i+1) = bd.Ptank_ox_v(i)*(bd.Vgas_ox(i)/bd.Vgas_ox(i+1))^gamma_He;
bd.Ptank_fuel_v(i+1) = bd.Ptank_fuel_v(i)*(bd.Vgas_fuel(i)/bd.Vgas_fuel(i+1))^gamma_He;

if cooling == "SI"
    [T_end_cooling,Delta_P_cool] = coolingIterative(bd.mdot_fuel_v(i), ...
                engineSize,engineSizeCEA,coolingSize);
    bd.cooling.T_end_cooling(i+1) = T_end_cooling;
    bd.cooling.Delta_P_cool(i+1) = Delta_P_cool;
    K_cool_fuel = (1-coolingSize.constants.beta^4)/(2*rho_fuel*(coolingSize.constants.C_D_cooling*coolingSize.constants.A_tube)^2);
else
    K_cool_fuel = 0;
end

% Solving system
funs = @(x) [-x(1)+bd.cstar_v(i)*(x(2)+x(3))/(engineSize.diameters.Dthroat^2*pi/4);
             -bd.Ptank_fuel_v(i+1)+x(1)+(K_inj_fuel+K_idr_fuel+K_dyn_fuel+K_cool_fuel)*x(2)^2;
             -bd.Ptank_ox_v(i+1)+x(1)+(K_inj_ox+K_idr_ox+K_dyn_ox)*x(3)^2];

sol = fsolve(funs,[bd.P_cc(i),bd.mdot_fuel_v(i),bd.mdot_ox_v(i)],optimoptions('fsolve','Display','none'));
bd.P_cc(i+1) = sol(1);
bd.mdot_fuel_v(i+1) = sol(2);
bd.mdot_ox_v(i+1) = sol(3);

% Total mass flow rate
bd.mdot_v(i+1) = bd.mdot_ox_v(i+1)+bd.mdot_fuel_v(i+1);

% O/F ratio
bd.O_F(i+1) = bd.mdot_ox_v(i+1)/bd.mdot_fuel_v(i+1);

% Thrust
bd.T(i+1) = (bd.mdot_v(i+1)*bd.cT_v(i)*bd.cstar_v(i)-bd.P_e(i)*engineSize.Areas.Aexit)*...
            engineSize.twoDlosses.lambda+bd.P_e(i)*engineSize.Areas.Aexit;

end

