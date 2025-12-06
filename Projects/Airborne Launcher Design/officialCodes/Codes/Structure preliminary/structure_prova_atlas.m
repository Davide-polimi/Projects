clc; clear; close all;

%% FIRST ATLAS V

diam_first = 3.81;
m_prop_stage_first = 284089;
rho_fu_first = 820;
rho_ox_first = 1140;
oxidizer_first = 'LOX';
fuel_first = 'RP-1';
of_mix_first = 2.7;
AR_tank_first = 1.4142;
T_first = 4152 * 1e+03;

mass_inert_first_real = 21054;

[masses_first, lengths_first] = structure_preliminary ('First', diam_first, m_prop_stage_first, rho_ox_first, rho_fu_first, of_mix_first, oxidizer_first, fuel_first, AR_tank_first, T_first)

mass_wiring_first = 1.43 * 32.46;

mass_avionics_first = 0.2*350;

inert_mass_first = masses_first.m_oxidizer_assembly + masses_first.m_fuel_assembly + ...
    masses_first.m_thrust_structure + masses_first.m_engine + mass_wiring_first + ...
    mass_avionics_first + masses_first.m_first_step_aft_skirt + ... 
    masses_first.m_intertank + masses_first.m_interstage

%inert_mass_first = 1.8 * inert_mass_first;

percentage_first = (inert_mass_first / mass_inert_first_real - 1)*100

% eps_s1 = inert_mass_first / (m_prop_stage_first + inert_mass_first)

%% SECOND ATLAS V

diam_second = 3.05;
m_prop_stage_second = 20830;
rho_fu_second = 71;
rho_ox_second = 1140;
oxidizer_second = 'LOX';
fuel_second = 'LH2';
of_mix_second = 5;
AR_tank_second = 1.4142;
T_second = 99200;

mass_inert_second_real = 2243;

% [masses_second, lengths_second] = structure_preliminary (diam_second, m_prop_stage_second, rho_ox_second, rho_fu_second, of_mix_second, oxidizer_second, fuel_second, AR_tank_second, T_second)

[masses_second, lengths_second] = structure_preliminary ('Second', diam_second, m_prop_stage_second, rho_ox_second, rho_fu_second, of_mix_second, oxidizer_second, fuel_second, AR_tank_second, T_second);

mass_wiring_second = 1.43 * 12.68;

mass_avionics_second = 0.8 * 350;

inert_mass_second = masses_second.m_oxidizer_assembly + masses_second.m_fuel_assembly + ...
    masses_second.m_engine + mass_wiring_second + ...
    mass_avionics_second + masses_second.m_intertank + ... 
    masses_second.m_second_step_aft_skirt + masses_second.m_forward_skirt

% inert_mass_second = 1.8 * inert_mass_second;

percentage_second = (inert_mass_second / mass_inert_second_real - 1) * 100

% eps_s2 = inert_mass_second / (m_prop_stage_second + inert_mass_second)

% masses_second.m_thrust_structure + masses_second.m_engine;
