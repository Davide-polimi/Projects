clearvars
close all
clc
%add_settings

% Preliminary geometry code
% Ref: slide 10 of "06 Mass, Forces, and Structures Part 1"

LsuD_LauncherOne = 12.38;

% DATA form LIQUIDLIQUID.m
rho_F = 810;    % FUEL density [kg/m^2]
rho_OX = 1142;  % OX density [kg/m^2]
M_F = [1.135354109843789e+04 5.141502913463922e+03];    % FUEL mass [kg]
M_OX = [4.205015221643661e+03 1.904260338319971e+03];   % OX mass [kg]
V_F = M_F / rho_F;      % FUEL volume [kg]
V_OX = M_OX / rho_OX;   % OX volume [kg]

% BILIQUID, 2 stages, tandem hemispherical tanks
D = [1.8 1.5];  % Diamter of stages [m]
V = [V_OX(1), V_F(1);
     V_OX(2), V_F(2)];

% DATA from ATLASV
D = [3.810000000000000, 3.050000000000000];
V = [2.257014191482911e+02, 1.374872655209796e+02; 15.226608187134506,  48.896713615023470];

AR = sqrt(2);

[Ltot,LsuD,h_stage_tot,L,h_stage,h_tank] = launcher_preliminary_lengths(D, V, AR);
Ltot
LsuD
