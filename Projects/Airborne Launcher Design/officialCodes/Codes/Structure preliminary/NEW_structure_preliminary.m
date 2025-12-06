function [MASSES, DIMENSIONS, PARAMETERS] = NEW_structure_preliminary (INPUT)

% NEW_structure_preliminary computes the bla bla
%
% INPUT:
% INPUT.n_stages            [double]    [scalar]
% INPUT.propulsion_type     [char]      [vector]
% INPUT.type_oxidizer       [char]      [vector]
% INPUT.type_fuel           [char]      [vector]
% INPUT.type_grain          [char]      [vector]
% INPUT.m_prop              [double]    [vector]
% INPUT.of_mix              [double]    [vector]
% INPUT.rho_fuel            [double]    [vector]
% INPUT.rho_ox              [double]    [vector]
% INPUT.rho_grain           [double]    [vector]
% INPUT.common_bulkhead     [double]    [vector, 0-1-2]
% INPUT.AR_tank             [double]    [vector]
% INPUT.T                   [double]    [vector]
% INPUT.eps_nozzle          [double]    [vector]
% INPUT.p_chamber           [double]    [vector]
% INPUT.m_payload           [double]    [scalar]
% INPUT.fairing             [struct]
%       .fairing.diameter       [double]    [scalar]
%       .fairing.length         [double]    [scalar]
%       .fairing.surface        [double]    [scalar]
%       .fairing.areal_density  [double]    [scalar]
%       .fairing.material       [char]      [scalar]

%% UNIFORMITY OF INPUT STRINGS
% In this section the INPUT fields that have type 'char' are uniformed to
% be used in the code

for kk = 1:INPUT.n_stages
    if ( strcmp (INPUT.propulsion_type(kk,:), 'Liquid') || ...
         strcmp (INPUT.propulsion_type(kk,:), 'liquid') || ...
         strcmp (INPUT.propulsion_type(kk,:), 'LRE') || ...
         strcmp (INPUT.propulsion_type(kk,:), 'lre') )
        INPUT.propulsion_type(kk,:) = 'Liquid';
    elseif ( strcmp (INPUT.propulsion_type(kk,:), 'Solid') || ...
             strcmp (INPUT.propulsion_type(kk,:), 'solid') || ...
             strcmp (INPUT.propulsion_type(kk,:), 'SRM') || ...
             strcmp (INPUT.propulsion_type(kk,:), 'srm') )
        INPUT.propulsion_type(kk,:) = 'Solid';
    end
end

for kk = 1:INPUT.n_stages
    if ( strcmp (INPUT.type_oxidizer(kk,:), 'LOX') || ...
         strcmp (INPUT.type_oxidizer(kk,:), 'lox') || ...
         strcmp (INPUT.type_oxidizer(kk,:), 'LO2') || ...
         strcmp (INPUT.type_oxidizer(kk,:), 'lo2') )
        INPUT.type_oxidizer(kk,:) = 'LOX';
    end
end

for kk = 1:INPUT.n_stages
    if ( strcmp (INPUT.type_fuel(kk,:), 'LH2') || ...
         strcmp (INPUT.type_fuel(kk,:), 'lh2') )
        INPUT.type_fuel(kk,:) = 'LH2';
    elseif ( strcmp (INPUT.type_fuel(kk,:), 'RP1') || ...
             strcmp (INPUT.type_fuel(kk,:), 'RP-1') || ...
             strcmp (INPUT.type_fuel(kk,:), 'rp1') || ...
             strcmp (INPUT.type_fuel(kk,:), 'rp-1') )
        INPUT.type_fuel(kk,:) = 'RP1';
    elseif ( strcmp (INPUT.type_fuel(kk,:), 'LCH4') || ...
             strcmp (INPUT.type_fuel(kk,:), 'CH4') || ...
             strcmp (INPUT.type_fuel(kk,:), 'lch4') || ...
             strcmp (INPUT.type_fuel(kk,:), 'ch4') )
        INPUT.type_fuel(kk,:) = 'LCH4';
    end
end

%% STRUCTURAL MASSES AND DIMENSIONS
% In this section the MASSES and DIMENSIONS of the various portions of the
% launcher are computed.
% Reference: Edberg and Costa

for kk = 1:INPUT.n_stages
    if strcmp(INPUT.propulsion_type(kk,:), 'Liquid')
        MASSES.fuel(kk) = INPUT.m_prop(kk) * (1 / (1 + INPUT.of_mix(kk)));
        MASSES.ox(kk) = INPUT.m_prop(kk) * (INPUT.of_mix(kk) / (1 + INPUT.of_mix(kk)));
    
        DIMENSIONS.vol_fuel(kk) = MASSES.fuel(kk) / INPUT.rho_fuel(kk);
        DIMENSIONS.vol_ox(kk) = MASSES.ox(kk) / INPUT.rho_ox(kk);

%% ---------- TANKS, INTERTANKS, THRUST (ENGINE AND STRUCTURE) ----------
% In this section, the MASSES and DIMENSIONS of the tanks and intertanks
% are computed. Furthermore, the MASSES of the engine and of the thrust 
% structure are computed.
        
        switch INPUT.common_bulkhead(kk)
            case 0
            % If common_bulkhead of the stage is 0, it means that the
            % oxidizer and fuel tanks are separated. In this case, there is
            % an additional length and mass dictated by the presenze of an
            % intertank structure. Please consider that the oxidizer tank
            % is the one that senses (in the code) the presence of the
            % common bulkhead option, since the fuel tank is always
            % considered as a usual tank with two domes.

            [DIMENSIONS.l_tank_ox(kk), DIMENSIONS.h_cyl_ox(kk), DIMENSIONS.h_dome_ox(kk), ...
             DIMENSIONS.A_tank_ox(kk), DIMENSIONS.vol_tank_ox(kk)] = ...
             tank_function (INPUT.diameter(kk), DIMENSIONS.vol_ox(kk), INPUT.AR_tank(kk), INPUT.common_bulkhead(kk));
            
            [DIMENSIONS.l_tank_fuel(kk), DIMENSIONS.h_cyl_fuel(kk), DIMENSIONS.h_dome_fuel(kk), ...
             DIMENSIONS.A_tank_fuel(kk), DIMENSIONS.vol_tank_fuel(kk)] = ...
             tank_function (INPUT.diameter(kk), DIMENSIONS.vol_fuel(kk), INPUT.AR_tank(kk), 0);
            
            DIMENSIONS.l_intertank(kk) = 1/4*INPUT.diameter(kk) + ...
                                         2*DIMENSIONS.h_dome_ox(kk);
            case 1
            % If common_bulkhead of the stage is 1, it means that the
            % OXIDIZER and fuel tanks share a dome. In this case, there is
            % a reduction in mass since there is not an intertank
            % structure, since the cylindrical part of the oxidizer tank is
            % in "direct contact" with the cylindrical part of the fuel
            % tank. Please consider that the OXIDIZER tank is the one that 
            % senses (in the code) the presence of the common bulkhead 
            % option.
            [DIMENSIONS.l_tank_ox(kk), DIMENSIONS.h_cyl_ox(kk), DIMENSIONS.h_dome_ox(kk), ...
             DIMENSIONS.A_tank_ox(kk), DIMENSIONS.vol_tank_ox(kk)] = ...
             tank_function (INPUT.diameter(kk), DIMENSIONS.vol_ox(kk), INPUT.AR_tank(kk), INPUT.common_bulkhead(kk));
            
            [DIMENSIONS.l_tank_fuel(kk), DIMENSIONS.h_cyl_fuel(kk), DIMENSIONS.h_dome_fuel(kk), ...
             DIMENSIONS.A_tank_fuel(kk), DIMENSIONS.vol_tank_fuel(kk)] = ...
             tank_function (INPUT.diameter(kk), DIMENSIONS.vol_fuel(kk), INPUT.AR_tank(kk), 0);
            
            DIMENSIONS.l_intertank(kk) = 0;

            case 2
            % If common_bulkhead of the stage is 2, it means that the
            % oxidizer and FUEL tanks share a dome. In this case, there is
            % a reduction in mass since there is not an intertank
            % structure, since the cylindrical part of the oxidizer tank is
            % in "direct contact" with the cylindrical part of the fuel
            % tank. Please consider that the FUEL tank is the one that 
            % senses (in the code) the presence of the common bulkhead 
            % option.
            [DIMENSIONS.l_tank_ox(kk), DIMENSIONS.h_cyl_ox(kk), DIMENSIONS.h_dome_ox(kk), ...
             DIMENSIONS.A_tank_ox(kk), DIMENSIONS.vol_tank_ox(kk)] = ...
             tank_function (INPUT.diameter(kk), DIMENSIONS.vol_ox(kk), INPUT.AR_tank(kk), 0);
            
            [DIMENSIONS.l_tank_fuel(kk), DIMENSIONS.h_cyl_fuel(kk), DIMENSIONS.h_dome_fuel(kk), ...
             DIMENSIONS.A_tank_fuel(kk), DIMENSIONS.vol_tank_fuel(kk)] = ...
             tank_function (INPUT.diameter(kk), DIMENSIONS.vol_fuel(kk), INPUT.AR_tank(kk), INPUT.common_bulkhead(kk));
            
            DIMENSIONS.l_intertank(kk) = 0;
                
        end

        if strcmp(INPUT.type_oxidizer(kk,:), 'LOX')
            % Please consider that the MER related to the mass of the tank
            % does not take into account some characteristics of the tank
            % (e.g. no thickness considered, no material considered, no 
            % shape considered)
            MASSES.tank_ox(kk) = 0.0107 * MASSES.ox(kk);
            MASSES.ins_ox(kk) = 1.123 * DIMENSIONS.A_tank_ox(kk);

        else
            MASSES.tank_ox(kk) = 12.16 * DIMENSIONS.vol_ox(kk);
            MASSES.ins_ox(kk) = 0;
        end

        if strcmp(INPUT.type_fuel(kk,:), 'LH2')
            % Please consider that the MER related to the mass of the tank
            % does not take into account some characteristics of the tank
            % (e.g. no thickness considered, no material considered, no 
            % shape considered)
            MASSES.tank_fuel(kk) = 0.128 * MASSES.fuel(kk);
            MASSES.ins_fuel(kk) = 2.88 * DIMENSIONS.A_tank_fuel(kk);
        
        elseif strcmp(INPUT.type_fuel(kk,:), 'RP1')
            % Please consider that the MER related to the mass of the tank
            % does not take into account some characteristics of the tank
            % (e.g. no thickness considered, no material considered, no 
            % shape considered)
            MASSES.tank_fuel(kk) = 0.0148 * MASSES.fuel(kk);
            MASSES.ins_fuel(kk) = 0;
        
        elseif strcmp(INPUT.type_fuel(kk,:), 'LCH4')
            MASSES.tank_fuel(kk) = 12.16 * DIMENSIONS.vol_fuel(kk);
            MASSES.ins_fuel(kk) = 0.98 * DIMENSIONS.A_tank_fuel(kk);

        else
            MASSES.tank_fuel(kk) = 12.16 * DIMENSIONS.vol_fuel(kk);
            MASSES.ins_fuel(kk) = 0;
        end

        MASSES.engine(kk) = (7.81*1e-4)*INPUT.T(kk) + ... 
                              (3.37*1e-5)*INPUT.T(kk)*sqrt(INPUT.eps_nozzle(kk)) ...
                              + 59;

        MASSES.thrust_structure(kk) = 2.55*1e-4 * INPUT.T(kk);

    elseif strcmp(INPUT.propulsion_type(kk,:), 'Solid')
        error ('SRM not implemented yet :(');
            %MASSES.SRM_case(kk) = 0.135 * INPUT.m_prop(kk);
    end
end

%% ---------- SKIRTS ----------
% In this section, the MASSES and DIMENSIONS of the skirts element are
% computed. Note that these components are not pressurized. 

DIMENSIONS.l_aft_skirt(1) = INPUT.diameter(1);

for kk = 2:INPUT.n_stages
    if strcmp (INPUT.propulsion_type(kk,:), 'Liquid')
        DIMENSIONS.l_aft_skirt(kk) = 1/3*INPUT.diameter(kk) + DIMENSIONS.h_dome_ox(kk);
    end
end

DIMENSIONS.l_forward_skirt = 1/3*INPUT.diameter(end) + DIMENSIONS.h_dome_fuel(end);
MASSES.forward_skirt = 13.3*pi * INPUT.diameter(end) * DIMENSIONS.l_forward_skirt;

%% ---------- INTERSTAGE ----------
% In this section the MASSES and DIMENSIONS of the interstage are computed.
% Please consider that the interstage is indexed as being part of the stage
% below it (so, the interstage between 1st and 2nd stage is indexed in
% position 1, the interstage between 2nd and 3rd stage is indexed in
% position 2, ecc.). Note that this component is not pressurized.

if INPUT.n_stages > 1
    for kk = 2:(INPUT.n_stages)
        if strcmp (INPUT.propulsion_type(kk,:), 'Liquid')
            switch INPUT.AR_tank(kk)
                case 1
                    DIMENSIONS.l_interstage(kk-1) = 5/4*INPUT.diameter(kk-1);
                    apotema = abs(0.5*INPUT.diameter(kk-1) - 0.5*INPUT.diameter(kk));
                    apotema = sqrt(apotema^2 + DIMENSIONS.l_interstage(kk-1)^2);
                    MASSES.interstage(kk-1) = ... 
                        13.3*pi * (0.5*INPUT.diameter(kk-1) + 0.5*INPUT.diameter(kk)) * apotema;
                otherwise
                    DIMENSIONS.l_interstage(kk-1) = INPUT.diameter(kk-1);
                    apotema = abs(0.5*INPUT.diameter(kk-1) - 0.5*INPUT.diameter(kk));
                    apotema = sqrt(apotema^2 + DIMENSIONS.l_interstage(kk-1)^2);
                    MASSES.interstage(kk-1) = ... 
                        13.3*pi * (0.5*INPUT.diameter(kk-1) + 0.5*INPUT.diameter(kk)) * apotema;
            end
        elseif strcmp (INPUT.propulsion_type(kk,:), 'Solid')
            error ('SRM not implemented yet :(\n');
        end
    end
end

MASSES.interstage(INPUT.n_stages) = 0;

for kk = 1:INPUT.n_stages
    MASSES.aft_skirt(kk) = 13.3*pi * INPUT.diameter(kk) * DIMENSIONS.l_aft_skirt(kk);
    MASSES.intertank(kk) = 13.3*pi * INPUT.diameter(kk) * DIMENSIONS.l_intertank(kk);
end

%% ---------- FAIRING ----------

if ~isfield(INPUT.fairing, 'diameter')
    DIMENSIONS.fairing.diameter = INPUT.diameter(end);
else
    DIMENSIONS.fairing.diameter = INPUT.fairing.diameter;
end

if ~isfield(INPUT.fairing, 'length')
    DIMENSIONS.fairing.length = 2*INPUT.diameter(end);
else
    DIMENSIONS.fairing.length = INPUT.fairing.length;
end

if ~isfield(INPUT.fairing, 'surface')
    R = DIMENSIONS.fairing.diameter / 2;
    h = DIMENSIONS.fairing.length;
    DIMENSIONS.fairing.surface = pi*R*sqrt(R^2 + h^2);
else
    DIMENSIONS.fairing.surface = INPUT.fairing.surface;
end

if ~isfield(INPUT.fairing, 'areal_density')
    if strcmp (INPUT.fairing.material, 'Metal')
        MASSES.fairing = 13.3 * DIMENSIONS.fairing.surface;
    elseif strcmp (INPUT.fairing.material, 'Composite')
        MASSES.fairing = 9.89 * DIMENSIONS.fairing.surface;
    end
else
    MASSES.fairing = INPUT.areal_density * DIMENSIONS.fairing.surface;
end

MASSES.PAF = 0.0755*INPUT.m_payload + 50;

%% ---------- SUM UP ----------

for kk = 1:INPUT.n_stages
    if strcmp (INPUT.propulsion_type(kk,:), 'Liquid')
        switch INPUT.common_bulkhead(kk)
            case 0
            DIMENSIONS.l_stage(kk) = DIMENSIONS.l_aft_skirt(kk) + ...
                                     DIMENSIONS.h_cyl_ox(kk) + ...
                                     DIMENSIONS.l_intertank(kk) + ...
                                     DIMENSIONS.h_cyl_fuel(kk);
            case 1
            DIMENSIONS.l_stage(kk) = DIMENSIONS.l_aft_skirt(kk) + ...
                                     DIMENSIONS.h_cyl_ox(kk) + ...
                                     DIMENSIONS.l_intertank(kk) + ...
                                     DIMENSIONS.h_cyl_fuel(kk);
             case 2
            DIMENSIONS.l_stage(kk) = DIMENSIONS.l_aft_skirt(kk) + ...
                                     DIMENSIONS.h_cyl_ox(kk) + ...
                                     DIMENSIONS.l_intertank(kk) + ...
                                     DIMENSIONS.h_cyl_fuel(kk);

        end
    elseif strcmp (INPUT.propulsion_type(kk,:), 'Solid')
        error ('SRM not implemented yet :(\n');
    end
end

DIMENSIONS.l_stage(end) = DIMENSIONS.l_stage(end) + DIMENSIONS.l_forward_skirt;

DIMENSIONS.l_total = sum(DIMENSIONS.l_stage) + sum(DIMENSIONS.l_interstage) + ...
                     DIMENSIONS.fairing.length;

MASSES.wiring = 1.43 * DIMENSIONS.l_total;

if INPUT.m_payload < 1000
    MASSES.avionics = 75;
else
    MASSES.avionics = 350;
end

for kk = 1:INPUT.n_stages
    if strcmp (INPUT.propulsion_type(kk,:), 'Liquid')
        MASSES.inert_stage(kk) = MASSES.thrust_structure(kk) + ...
                                 MASSES.engine(kk) + ...
                                 MASSES.tank_ox(kk) + ...
                                 MASSES.ins_ox(kk) + ...
                                 MASSES.tank_fuel(kk) + ...
                                 MASSES.ins_fuel(kk) + ...
                                 MASSES.aft_skirt(kk) + ...
                                 MASSES.intertank(kk) + ...
                                 MASSES.wiring*(DIMENSIONS.l_stage(kk)/DIMENSIONS.l_total);
    end
    if ~isequal(kk, INPUT.n_stages)
        MASSES.inert_stage(kk) = MASSES.inert_stage(kk) + ...
                                 MASSES.avionics * (0.2/INPUT.n_stages);
    else 
        MASSES.inert_stage(kk) = MASSES.inert_stage(kk) + ...
                                 MASSES.avionics * 0.8;
    end
end

for kk = 1:INPUT.n_stages
    MASSES.structural(kk) = MASSES.inert_stage(kk) + ...
                            MASSES.interstage(kk);
end

for kk = 1:INPUT.n_stages
    PARAMETERS.eps_str(kk) = MASSES.structural(kk) / (MASSES.structural(kk) + ...
                             INPUT.m_prop(kk));
end

end