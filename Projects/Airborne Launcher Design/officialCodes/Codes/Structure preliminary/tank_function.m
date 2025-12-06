function [l_tank, h_cyl, h_dome, A_tank, V_tank] = tank_function (diam, vol_prop, AR, bulkhead_option, sphere_option)

% tank_function provides the dimensions of a tank. By default, a cylindric + dome
% shape is computed. Further option for spherical tank can be specified.
%
% INPUT:
% - diam            [double]    [m]     Diameter
% - vol_prop        [double]    [m3]    Propellant volume
% - AR              [double]    [-]     Aspect ratio of the tank
% - bulkhead_option [double]    [scalar, 0-1]   Option for common bulkhead
% - sphere_option   [char]              Option to compute a spherical tank
%
% OUTPUT:
% - l_tank          [m]     Length of the tank
% - h_cyl           [m]     Length of the cylindrical part
% - h_dome          [m]     Length of the dome part
% - A_tank          [m2]    Area of the tank
% - V_tank          [m3]    Volume of the tank (CHECK: V_tank == vol_prop)

if nargin == 4
    R = diam / 2;

    switch bulkhead_option
        case 0
    
        h_cyl = 4 * vol_prop / (pi * diam^2) - 2*diam / (3 * AR);
        h_dome = R / AR;
        
        l_tank = h_cyl + 2*h_dome;
        
        E = sqrt(1 - 1/(AR^2));
        
        switch AR
            case 1 
                A_tank = 2*pi*R*h_cyl + 4*pi*R^2;
                V_tank = pi*R^2*h_cyl + 4/3*pi*R^3;
                    
            otherwise
                A_tank = 2*pi*R*h_cyl + 2*pi*R^2 * (1 + 1 / (2*E*AR^2) * log((1+E)/(1-E)));
                V_tank = pi*R^2*h_cyl + 4/3*pi*R^3/AR;
        
        end

        case 1

        h_cyl = 4 * vol_prop / (pi * diam^2) - diam / (3 * AR);
        h_dome = R / AR;
        
        l_tank = h_cyl + h_dome;
        
        E = sqrt(1 - 1/(AR^2));
        
        switch AR
            case 1 
                A_tank = 2*pi*R*h_cyl + 2*pi*R^2;
                V_tank = pi*R^2*h_cyl + 2/3*pi*R^3;
                    
            otherwise
                A_tank = 2*pi*R*h_cyl + pi*R^2 * (1 + 1 / (2*E*AR^2) * log((1+E)/(1-E)));
                V_tank = pi*R^2*h_cyl + 2/3*pi*R^3/AR;
        end

        case 2

        h_cyl = 4 * vol_prop / (pi * diam^2) - diam / (3 * AR);
        h_dome = R / AR;
        
        l_tank = h_cyl + h_dome;
        
        E = sqrt(1 - 1/(AR^2));
        
        switch AR
            case 1 
                A_tank = 2*pi*R*h_cyl + 2*pi*R^2;
                V_tank = pi*R^2*h_cyl + 2/3*pi*R^3;
                    
            otherwise
                A_tank = 2*pi*R*h_cyl + pi*R^2 * (1 + 1 / (2*E*AR^2) * log((1+E)/(1-E)));
                V_tank = pi*R^2*h_cyl + 2/3*pi*R^3/AR;
        
        end

    end

elseif nargin == 5
    if strcmp (sphere_option, 'sphere')
        R = diam/2;
        h_cyl = 0;
        h_dome = R;
        l_tank = h_cyl + 2*h_dome;
        A_tank = 4*pi*R^2;
        V_tank = 4/3*pi*R^3;
    else
        R = diam / 2;
    
    h_cyl = 4 * vol_prop / (pi * diam^2) - 2*diam / (3 * AR);
    h_dome = R / AR;
    
    l_tank = h_cyl + 2*h_dome;
    
    E = sqrt(1 - 1/(AR^2));
    
    switch AR
        case 1 
            A_tank = 2*pi*R*h_cyl + 4*pi*R^2;
            V_tank = pi*R^2*h_cyl + 4/3*pi*R^3;
             
        otherwise
            A_tank = 2*pi*R*h_cyl + 2*pi*R^2 * (1 + 1 / (2*E*AR^2) * log((1+E)/(1-E)));
            V_tank = pi*R^2*h_cyl + 4/3*pi*R^3/AR;
    
    end

    end

end

end