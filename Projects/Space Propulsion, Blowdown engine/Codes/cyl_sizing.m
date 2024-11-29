function [tank,Gas] = cyl_sizing(B_opt,tank,r_avail,h_avail,engineSize,material, ...
                                Gas,rho_ox,rho_fuel,temp_LOX,temp_RP1)

% Nota: pipe_empty è lo spazio libero per il passaggio delle pipes
% Nota: in = interno
% Nota: r_ox_in =       x(1)
% Nota: th_tank_ox =    x(2)
% Nota: h_tank_ox =     x(3)
% Nota: r_fuel_in =     x(4)
% Nota: th_tank_fuel =  x(5)
% Nota: h_tank_fuel =   x(6)

L_conv = engineSize.lengths.Lconv;
L_c = engineSize.lengths.Lchamber;
r_c = engineSize.diameters.Dchamber/2;
r_t = engineSize.diameters.Dthroat/2;

sw_ox = material.sigma_y.crio/1.33; % Nota: 1.33 arriva da pag. 335 del design of liquid propellant engines
                                    % 1971 (huzel), 19710019929.pdf nel caso servisse
sw_fuel = material.sigma_y.amb/1.33;

sizing = @(x) [-r_avail+x(1)+x(2); 
               -r_avail+tank.pipes_empty+x(4)+x(5);
               -x(3)-2*x(2)+(0.8*r_avail^2*h_avail*pi-L_c*pi*r_c^2-pi/3*L_conv*((r_c+r_t)/2)^2-(x(4)+x(5))^2*pi*(x(6)+2*x(5)))/((x(1)+x(2))^2*pi);
               -x(6)-2*x(5)+(x(3)+2*x(2))*(x(1)+x(2))^2/(tank.V_ratio*(x(4)+x(5))^2);
               -x(2)+tank.P_ox*x(1)/sw_ox;
               -x(5)+tank.P_fuel*x(4)/sw_fuel];

% Nota: first guess
% 0.48 = raggio ox
% 0.007 th ox
% 0.8 altezza ox

x0 = [0.48, 0.007, 0.8, 0.4, 0.009, 0.8];

sol = fsolve(sizing, x0, optimoptions("fsolve", "Display","none"));

tank.r_ox_in = sol(1);
tank.th_ox = sol(2);
tank.h_ox = sol(3);
tank.r_fuel_in = sol(4);
tank.th_fuel = sol(5);
tank.h_fuel = sol(6);

tank.V_ox_tot = tank.r_ox_in^2*pi*tank.h_ox;        % Volume INTERNO oxidizer
tank.V_fuel_tot = tank.r_fuel_in^2*pi*tank.h_fuel;  % Volume INTERNO fuel


tank.Volume_check = ((tank.r_ox_in+tank.th_ox)^2*pi*(tank.h_ox+2*tank.th_ox)+ ...
                    (tank.r_fuel_in+tank.th_fuel)^2*pi*(tank.h_fuel+2*tank.th_fuel)+ ...
                    L_c*pi*r_c^2+pi/3*L_conv*((r_c+r_t)/2)^2)/(r_avail^2*h_avail*pi) ;

tank.w_ox = material.rho.crio*((tank.r_ox_in+tank.th_ox)^2*pi*(tank.h_ox+2*tank.th_ox)- ...
            tank.r_ox_in^2*pi*tank.h_ox);
tank.w_fuel = material.rho.amb*((tank.r_fuel_in+tank.th_fuel)^2*pi*(tank.h_fuel+2*tank.th_fuel)- ...
            tank.r_fuel_in^2*pi*tank.h_fuel);

% Blow-down optimization 
tank.B = B_opt;

tank.Vgas_ox_fin = tank.V_ox_tot;
tank.Vgas_ox_in = tank.Vgas_ox_fin/(tank.B^(1/Gas.k));
Gas.V_ox_in = tank.Vgas_ox_in;

tank.V_ox_in = tank.Vgas_ox_in*(tank.B^(1/Gas.k)-1);           % Volume oxidizer

tank.m_ox_in = tank.V_ox_in*rho_ox;
tank.mgas_ox_in = tank.Vgas_ox_in*tank.P_max/(Gas.R*temp_LOX);
Gas.mgas_ox = tank.mgas_ox_in;

tank.Vgas_fuel_fin = tank.V_fuel_tot;
tank.Vgas_fuel_in = tank.Vgas_fuel_fin/(tank.B^(1/Gas.k));
Gas.V_fuel_in = tank.Vgas_fuel_in;

tank.V_fuel_in = tank.Vgas_fuel_in*(tank.B^(1/Gas.k)-1);       % Volume oxidizer

tank.m_fuel_in = tank.V_fuel_in*rho_fuel;
tank.mgas_fuel_in = tank.Vgas_fuel_in*tank.P_max/(Gas.R*temp_RP1);
Gas.mgas_fuel = tank.mgas_fuel_in;

end