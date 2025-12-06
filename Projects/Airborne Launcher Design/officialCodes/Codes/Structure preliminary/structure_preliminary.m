function [masses, lengths] = structure_preliminary (which_stage, diam, m_prop_stage, rho_ox, rho_fu, of_mix, oxidizer, fuel, AR_tank, T)

% structure_preliminary provides the structs masses and lengths by a
% preliminary estimation of some components masses.
% The relations are the MERs reported by Edberg and Costa book in chapter
% 7.
% Future improvements: refinement of the single components and other
% components not yet considered to be added. 
%
% INPUT:
% - which_stage     [char]  Specify the stage of interest
% - diam            [m]     Diameter of the stage
% - m_prop_stage    [kg]    Propellant mass of the stage
% - rho_ox          [kg/m3] Oxidizer density
% - rho_fu          [kg/m3] Fuel density
% - of_mix          [-]     O/F ratio
% - oxidizer        [char]  Oxidizer name (e.g. 'LOX')
% - fuel            [char]  Fuel name (e.g. 'LH2', 'RP-1')
% - AR_tank         [-]     Aspect Ratio of the tank
% - T               [N]     Thrust
%
% OUTPUT:
% - masses  [struct]    Struct of the masses
% - lengths [struct]    Struct of the dimensions (lengths, areas, volumes)

masses.m_fu = m_prop_stage * (1 / (1 + of_mix));
masses.m_ox = m_prop_stage * (of_mix / (1 + of_mix));

lengths.vol_fu = masses.m_fu / rho_fu;
lengths.vol_ox = masses.m_ox / rho_ox;

% ---------- TANK DESIGN & MASS ----------

[lengths.l_tank_ox, lengths.h_cyl_ox, lengths.h_dome_ox, lengths.A_tank_ox, lengths.vol_tank_ox] = tank_function (diam, lengths.vol_ox, AR_tank);
[lengths.l_tank_fu, lengths.h_cyl_fu, lengths.h_dome_fu, lengths.A_tank_fu, lengths.vol_tank_fu] = tank_function (diam, lengths.vol_fu, AR_tank);

switch oxidizer
    case {'LOX', 'L02'}
        masses.m_ox_tank = 0.0107 * masses.m_ox;
        masses.m_ox_ins = 1.123 * lengths.A_tank_ox;

    otherwise
        error ('Oxidizers different from LOX (LO2) not implemented yet :(');
end

switch fuel
    case 'LH2'
        masses.m_fu_tank = 0.128 * masses.m_fu;
        masses.m_fu_ins = 2.88 * lengths.A_tank_fu;

    case {'RP1', 'RP-1'}
        masses.m_fu_tank = 0.0148 * masses.m_fu;
        masses.m_fu_ins = 0;

    otherwise
        error ('Fuels different from LH2 and RP-1 not implemented yet :(');
end

masses.m_oxidizer_assembly = masses.m_ox_tank + masses.m_ox_ins;
masses.m_fuel_assembly = masses.m_fu_tank + masses.m_fu_ins;

% ---------- ENGINE MASS ----------

T_lbs = units_conversion(T, 'newton_to_lbs');

if strcmp (oxidizer, 'LOX') && (strcmp (fuel, 'RP-1') || strcmp (fuel, 'RP1'))
    m_engine_lbs = 1200 + 0.00003 * (T_lbs^1.4);
    m_engine_kg = units_conversion(m_engine_lbs, 'lbs_to_kg');
elseif strcmp (oxidizer, 'LOX') && strcmp (fuel, 'LH2')
    m_engine_lbs = 150 + 0.086 * (T_lbs^0.86);
    m_engine_kg = units_conversion(m_engine_lbs, 'lbs_to_kg');
end

masses.m_engine = m_engine_kg;

% masses.m_engine = (7.81*1e-4)*T + (3.37*1e-5)*T*sqrt(eps_noz) + 59;

% ---------- THRUST STRUCTURE MASS ----------
switch which_stage
    case {'First', 'first', '1'}
        masses.m_thrust_structure = 2.55 * (1e-04) * T;
    otherwise
        masses.m_thrust_structure = [];
        disp ('No thrust structure for the selected stage\n');
end

% ---------- UNPRESSURIZED PARTS MASS ----------

% Considere that, to compute the areas of the single components (since the 
% MER for this part is generically 13.3 kg/m^2), the diameter is one of the
% function inputs and the length is retrieved from the general
% relationships of the missile length (Edberg and Costa pag. 479)

masses.m_first_step_aft_skirt = 13.3 * (pi*diam) * diam;
masses.m_intertank = 13.3 * (pi*diam) * (1/4*diam + lengths.h_dome_ox + lengths.h_dome_fu);

switch AR_tank
    case 1
        masses.m_interstage = 13.3 * (pi*diam) * (5/4*diam);
    otherwise
        masses.m_interstage = 13.3 * (pi*diam) * diam;
end

masses.m_second_step_aft_skirt = 13.3 * (pi*diam) * (1/3*diam + lengths.h_dome_ox);
masses.m_forward_skirt = 13.3 * (pi*diam) * (1/3*diam + lengths.h_dome_ox);

switch which_stage
    case {'First', 'first', '1'}
        masses.m_second_step_aft_skirt = [];
        masses.m_forward_skirt = [];
    case {'Second', 'second', '2'}
        masses.m_first_step_aft_skirt = [];
        masses.interstage = [];
end