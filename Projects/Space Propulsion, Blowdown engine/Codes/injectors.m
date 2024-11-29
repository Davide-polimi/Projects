function [geom_inj] = injectors(dP,mdot_ox,mdot_fu,rho_fu,rho_ox,Cd_ox,Cd_fu)

% INJECTORS SIZING

geom_inj.N_ox = 3;                  % 3 triplets
geom_inj.N_fu = 2*geom_inj.N_ox;    % triplet

geom_inj.A_fu_tot = mdot_fu/(Cd_fu*sqrt(2*rho_fu*dP));
geom_inj.A_ox_tot = mdot_ox/(Cd_ox*sqrt(2*rho_ox*dP));

geom_inj.A_fu_single = geom_inj.A_fu_tot/geom_inj.N_fu;
geom_inj.D_fu_single = 2*sqrt(geom_inj.A_fu_single/pi);

geom_inj.A_ox_single = geom_inj.A_ox_tot/geom_inj.N_ox; 
geom_inj.D_ox_single = 2*sqrt(geom_inj.A_ox_single/pi);

% angle between oxidizer jets
geom_inj.beta_imp = pi/3;

% impingement point (constraint on injection diameter)
geom_inj.L_imp = geom_inj.D_fu_single*5;              

% check on velocity
geom_inj.vel_ox = mdot_ox/(geom_inj.N_ox*geom_inj.A_ox_single*rho_ox);
geom_inj.vel_fu = mdot_fu/(geom_inj.N_fu*geom_inj.A_fu_single*rho_fu);

