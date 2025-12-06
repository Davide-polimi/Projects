clear; close; clc;
% -------------PARAMETRO DA VARIARE: DIAMETRO DEL SATELLITE----------------
% Ovviamente cambiando il diametro di base del satellite cambia il diametro del fairing
% va settato un valore plausibile (1 : 1.5 m circa)
diameterSat = 1.3; 

% dati:
density = 200; % s/c average density (0.2 g/cm^3)
mass = 400;    % Max Payload Mass [kg]

volumeSat = mass/density;
volumeIntFairing = 1.5*volumeSat; % DA CAPIRE OVERDIMENSIONAMENTO FAIRING
lengthSat = volumeSat/(pi*diameterSat^2 * 0.25);

%% LINEAR REGRESSION:
% Dati da lanciatori esistenti:
L_sat = [1.11; 0.717; 1.201; 2.28; 2.88; 4.888];
L_fair = [2.138; 1.9246; 2.230; 3.63; 2.88+1.524+0.876; 10.184];
D_fair = [1.152 + 0.2; 1.2036; 1.194 + 0.2; 1.36+0.2; 2.511+0.2; 4.572];


c1 = polyfit(L_sat, L_fair, 1);
L_sat_vec = linspace(0, 5, 1000);
L_fair_vec1 = polyval(c1, L_sat_vec);

c2 = polyfit(L_fair, D_fair, 1);
L_fair_vec2 = linspace(0, 11, 1000);
D_fair_vec = polyval(c2, L_fair_vec2);


figure(1)
plot(L_sat, L_fair, 'x', 'LineWidth', 1.5)
hold on
plot(L_sat_vec, L_fair_vec1, 'LineWidth', 1.5)
xlabel('Satellite Length [m]')
ylabel('Fairing Length [m]')


figure(2)
plot(L_fair, D_fair, 'x', 'LineWidth', 1.5)
hold on
plot(L_fair_vec2, D_fair_vec, 'LineWidth', 1.5)
xlabel('Fairing Length [m]')
ylabel('Fairing Diameter [m]')


%%
ourLsat = lengthSat;
ourLfairing = polyval(c1, ourLsat);
ourDfairing = polyval(c2, ourLfairing);

fprintf('Fairing Length: %.2f [m]\n', ourLfairing);
fprintf('Fairing Diameter: %.2f [m]\n', ourDfairing);
