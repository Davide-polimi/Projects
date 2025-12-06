%% from Payload mass to fairing diameter
clear 
close all
clc
% Dati noti
Mp = [7392; 285; 443; 200; 2400; 1140]; %ATLAS V 401 LEO, LAUNCHER ONE, PEGASUS XL, FALCON 1, VEGA C, VEGA
Dfair = [4.2; 1.37; 1.18; 1.524; 3; 2.6]; 

% Grado del polinomio
grado = 1; % ad esempio, retta (grado 1)

% Fit dei minimi quadrati
p = polyfit(Mp, Dfair, grado);

% Calcola il valore interpolato (predetto)
Mpq = 400; % punto di query
Dfairq = polyval(p, Mpq); % valore predetto dal modello

% Visualizzazione
xp = linspace(min(Mp), max(Mp), 100); % punti per grafico
yp = polyval(p, xp); % valutazione del polinomio

plot(Mp, Dfair, 'o', xp, yp, '-'); % grafico dati e interpolazione
xlabel('Mpayload'); ylabel('diameter');
title(['Interpolazione di grado ', num2str(grado)]);
legend('Dati', 'Fit polinomiale');

%% from Payload mass to fairing diameter
clear
close all
clc
% Dati noti
Mp = [7392; 200; 2400; 1140]; %ATLAS V 401 LEO, LAUNCHER ONE, PEGASUS XL, FALCON 1, VEGA C, VEGA
Mfair = [2305; 180; 860; 540]; % Fairing mass

% Grado del polinomio
grado = 1; % ad esempio, retta (grado 1)

% Fit dei minimi quadrati
p = polyfit(Mp, Mfair, grado);

% Calcola il valore interpolato (predetto)
Mpq = 250; % punto di query
Mfairq = polyval(p, Mpq); % valore predetto dal modello

% Visualizzazione
xp = linspace(min(Mp), max(Mp), 100); % punti per grafico
yp = polyval(p, xp); % valutazione del polinomio

plot(Mp, Mfair, 'o', xp, yp, '-'); % grafico dati e interpolazione
xlabel('Mpayload'); ylabel('Mfairing');
title(['Interpolazione di grado ', num2str(grado)]);
legend('Dati', 'Fit polinomiale');

%% from fairing diameter to fairing length
% Dati noti
Dfair = [4.2; 1.37; 1.18; 1.524; 3; 2.6]; %ATLAS V 401 LEO, LAUNCHER ONE, PEGASUS XL, FALCON 1, VEGA C, VEGA
Lfair = [12.68; 3.63; 2.13; 4.4; 9; 7.88]; 

% Grado del polinomio
grado = 1; % ad esempio, retta (grado 1)

% Fit dei minimi quadrati
p = polyfit(Dfair, Lfair, grado);

% Calcola il valore interpolato (predetto)
Lfairq = polyval(p, Dfairq); % valore predetto dal modello

% Visualizzazione
xp = linspace(min(Dfair), max(Dfair), 100); % punti per grafico
yp = polyval(p, xp); % valutazione del polinomio

plot(Dfair, Lfair, 'o', xp, yp, '-'); % grafico dati e interpolazione
xlabel('Fairing diameter'); ylabel('Length diameter');
title(['Interpolazione di grado ', num2str(grado)]);
legend('Dati', 'Fit polinomiale');


%% from payload mass to stages' thrust
%%%%%%1ST STAGE%%%%%%%
clear
close all 
clc
% Dati noti
Mp = [285; 443]; %LAUNCHER ONE, PEGASUS XL
T1stage = [327; 726]*1e3; 

% Grado del polinomio
grado = 1; % ad esempio, retta (grado 1)

% Fit dei minimi quadrati
p = polyfit(Mp, T1stage, grado);

% Calcola il valore interpolato (predetto)
Mpq = 250; % punto di query
T1stageq = polyval(p, Mpq); % valore predetto dal modello

% Visualizzazione
xp = linspace(min(Mp), max(Mp), 100); % punti per grafico
yp = polyval(p, xp); % valutazione del polinomio

plot(Mp, T1stage, 'o', xp, yp, '-'); % grafico dati e interpolazione
xlabel('Mpayload'); ylabel('Thrust 1st stage');
title(['Interpolazione di grado ', num2str(grado)]);
legend('Dati', 'Fit polinomiale');




%%%%%%2ND STAGE%%%%%%%

% Dati noti
Mp = [300; 285; 443; 200]; %Electron, LAUNCHER ONE, PEGASUS XL, FALCON 1
T2stage = [25.8; 22.241; 196; 30.7]*1e3; 

% Grado del polinomio
grado = 1; % ad esempio, retta (grado 1)

% Fit dei minimi quadrati
p = polyfit(Mp, T2stage, grado);

% Calcola il valore interpolato (predetto)
Mpq = 250; % punto di query
T2stageq = polyval(p, Mpq); % valore predetto dal modello

% Visualizzazione
xp = linspace(min(Mp), max(Mp), 100); % punti per grafico
yp = polyval(p, xp); % valutazione del polinomio

figure
plot(Mp, T2stage, 'o', xp, yp, '-'); % grafico dati e interpolazione
xlabel('Mpayload'); ylabel('Thrust 2st stage');
title(['Interpolazione di grado ', num2str(grado)]);
legend('Dati', 'Fit polinomiale');

%%%%%%3RD STAGE%%%%%%%

% Dati noti
Mp = [7392; 285; 443; 200]; %ATLAS V 401 LEO, LAUNCHER ONE, PEGASUS XL, FALCON 1, VEGA C, VEGA
T2stage = [198.4; 22.241; 196; 30.7]*1e3; 
Weight3rdstage = [];

% Grado del polinomio
grado = 1; % ad esempio, retta (grado 1)

% Fit dei minimi quadrati
p = polyfit(Mp, T2stage, grado);

% Calcola il valore interpolato (predetto)
Mpq = 250; % punto di query
T2stageq = polyval(p, Mpq); % valore predetto dal modello

% Visualizzazione
xp = linspace(min(Mp), max(Mp), 100); % punti per grafico
yp = polyval(p, xp); % valutazione del polinomio

figure
plot(Mp, T2stage, 'o', xp, yp, '-'); % grafico dati e interpolazione
xlabel('Mpayload'); ylabel('Thrust 2st stage');
title(['Interpolazione di grado ', num2str(grado)]);
legend('Dati', 'Fit polinomiale');


%% from payload mass to length over diameter

%L_total with D first stage

% Dati noti
% Mp = [7392; 285; 443; 200; 22800];
L = [21.336; 17.6; 21; 18; 70; 110.6; 30; 35; 58.3; 62];  %LAUNCHER ONE, PEGASUS XL, FALCON 1, ELECTRON,  FALCON 9 v1.2 
D = [1.8; 1.28; 1.68; 1.2;3.7; 10; 3; 3.4; 3.81; 5.4]; % First stage diameter, assumed constant


% Grado del polinomio
grado = 1; % ad esempio, retta (grado 1)

% Fit dei minimi quadrati
p = polyfit(D, L, grado);

% Calcola il valore interpolato (predetto)
D_guess = 1.7; % punto di query
L_guess = polyval(p, D_guess); % valore predetto dal modello

% Visualizzazione
xp = linspace(min(D), max(D), 100); % punti per grafico
yp = polyval(p, xp); % valutazione del polinomio

figure
plot(D, L, 'o', xp, yp, '-'); % grafico dati e interpolazione
xlabel('Mpayload'); ylabel('L/D 1st stage');
title(['Interpolazione di grado ', num2str(grado)]);
legend('Dati', 'Fit polinomiale');


%% from thrust to tb
%%%%%%1ST STAGE%%%%%%%

% Dati noti
%LAUNCHER ONE, PEGASUS XL
T1stage = [327; 726]*1e3; 
tb1 = [180; 68.6];

% Grado del polinomio
grado = 1; % ad esempio, retta (grado 1)

% Fit dei minimi quadrati
p = polyfit(T1stage, tb1, grado);

% Calcola il valore interpolato (predetto)
tb1q = polyval(p, T1stageq); % valore predetto dal modello

% Visualizzazione
xp = linspace(min(T1stage), max(T1stage), 100); % punti per grafico
yp = polyval(p, xp); % valutazione del polinomio

plot(T1stage, tb1, 'o', xp, yp, '-'); % grafico dati e interpolazione
xlabel('T1stage'); ylabel('tb');
title(['Interpolazione di grado ', num2str(grado)]);
legend('Dati', 'Fit polinomiale');




%%%%%%2ND STAGE%%%%%%%

% Dati noti
%LAUNCHER ONE, PEGASUS XL, FALCON 1, Electron
T2stage = [22.241; 196; 30.7; 25.8]*1e3;
tb2 = [360; 69.4; 378; 373];

% Grado del polinomio
grado = 1; % ad esempio, retta (grado 1)

% Fit dei minimi quadrati
p = polyfit(T2stage, tb2, grado);

% Calcola il valore interpolato (predetto)
tb2q = polyval(p, T2stageq); % valore predetto dal modello

% Visualizzazione
xp = linspace(min(T2stage), max(T2stage), 100); % punti per grafico
yp = polyval(p, xp); % valutazione del polinomio

figure
plot(T2stage, tb2, 'o', xp, yp, '-'); % grafico dati e interpolazione
xlabel('T2stage'); ylabel('tb');
title(['Interpolazione di grado ', num2str(grado)]);
legend('Dati', 'Fit polinomiale');





%% from payload mass to propellant mass and to structural mass index
%%%%%%1ST STAGE%%%%%%%
clear all
close all 
clc
% Dati noti
Mp = [285; 443]; %LAUNCHER ONE, PEGASUS XL
T1stage = [327; 726]*1e3; 

% Grado del polinomio
grado = 1; % ad esempio, retta (grado 1)

% Fit dei minimi quadrati
p = polyfit(Mp, T1stage, grado);

% Calcola il valore interpolato (predetto)
Mpq = 250; % punto di query
T1stageq = polyval(p, Mpq); % valore predetto dal modello

% Visualizzazione
xp = linspace(min(Mp), max(Mp), 100); % punti per grafico
yp = polyval(p, xp); % valutazione del polinomio

plot(Mp, T1stage, 'o', xp, yp, '-'); % grafico dati e interpolazione
xlabel('Mpayload'); ylabel('Thrust 1st stage');
title(['Interpolazione di grado ', num2str(grado)]);
legend('Dati', 'Fit polinomiale');



%%%%%%2nd STAGE%%%%%%%
% Dati noti
Mp = [285; 443]; %LAUNCHER ONE, PEGASUS XL
T1stage = [327; 726]*1e3; 

% Grado del polinomio
grado = 1; % ad esempio, retta (grado 1)

% Fit dei minimi quadrati
p = polyfit(Mp, T1stage, grado);

% Calcola il valore interpolato (predetto)
Mpq = 250; % punto di query
T1stageq = polyval(p, Mpq); % valore predetto dal modello

% Visualizzazione
xp = linspace(min(Mp), max(Mp), 100); % punti per grafico
yp = polyval(p, xp); % valutazione del polinomio

plot(Mp, T1stage, 'o', xp, yp, '-'); % grafico dati e interpolazione
xlabel('Mpayload'); ylabel('Thrust 1st stage');
title(['Interpolazione di grado ', num2str(grado)]);
legend('Dati', 'Fit polinomiale');