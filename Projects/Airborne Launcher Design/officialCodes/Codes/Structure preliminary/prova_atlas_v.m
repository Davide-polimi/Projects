clc; clear; close all;

ATLAS_V.n_stages = 2;
ATLAS_V.propulsion_type = ['Liquid'; 'Liquid'];
ATLAS_V.type_oxidizer = ['LOX'; 'LOX'];
ATLAS_V.rho_ox = [1140; 1140];
ATLAS_V.type_fuel = ['RP1'; 'LH2'];
ATLAS_V.rho_fuel = [820; 71];
ATLAS_V.of_mix = [2.72; 5.88];
ATLAS_V.m_prop = [284089; 20830];
ATLAS_V.diameter = [3.81; 3.05];
ATLAS_V.AR_tank = [1.4142; 1.4142];
ATLAS_V.T = [4152 * 1e+03; 99200];
ATLAS_V.eps_nozzle = [36.87; 100];
ATLAS_V.common_bulkhead = [0; 0];
ATLAS_V.m_payload = 4000;
ATLAS_V.fairing.material = 'Metal';
ATLAS_V.fairing.diameter = 4.2;
ATLAS_V.fairing.length = 12;

[MASSES, DIMENSIONS, PARAMETERS] = NEW_structure_preliminary(ATLAS_V);

ATLAS_V.inert_masses_real = [21054; 2243];
ATLAS_V.l_total_real = 32.46 + 2.52 + 0.65 + 12 + 12.68;

error_m_first = (MASSES.inert_stage(1) / ATLAS_V.inert_masses_real(1)) - 1;
error_m_second = (MASSES.inert_stage(2) / ATLAS_V.inert_masses_real(2)) - 1;
error_l_total = (DIMENSIONS.l_total / ATLAS_V.l_total_real) - 1;