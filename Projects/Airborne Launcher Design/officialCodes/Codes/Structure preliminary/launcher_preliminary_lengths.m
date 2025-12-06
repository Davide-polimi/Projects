function [Ltot,LsuD,h_stage_tot,L,h_stage,h_tank] = launcher_preliminary_lengths(D, V, AR)
% Preliminary geometry function
% Ref: slide 10,14 of "06 Mass, Forces, and Structures Part 1"
%
% INPUT:
% D         diameter vector D = [Dstage1, Dstage2, etc]                      [m]
% V         propellant volumes matrix: each row is a different stage, 
%           only two columns (one column for the fuel one for the ox)        [m^3]
%           example: V = [VoxStage1, VfStage1;
%                         VoxStage2, VfStage2]
% AR        tanks aspect ratio: 1 for hemispherical domes, other values
%           (sqrt(2) most used) for helliptical domes                        [-]
%
% OUTPUT:
% Ltot      total laucher length                                             [m]
% LsuD      Length/Diamter of the launcher                                   [-]
% h_stage   stages lengths matrix: each row composed of 
%           h_stage(i) = [h_aftSkirt h_cyl h_interTanks h_cyl h_interStage]  [m]
% h_tank    length of the tanks matrix (composed as the V matrix)            [m]
%
% NOTE:
%           Works only for biliquid, with turbopumps

n = length(D);
h_tank = zeros(n,2);
h_cyl = h_tank;
h_dome = h_tank;

for ii = 1:n
    for jj = 1:2
        [h_tank(ii,jj), h_cyl(ii,jj), h_dome(ii,jj)] = tank_function (D(ii),V(ii,jj), AR);
    end
end

h_interTanks = zeros(1,n);
h_aftSkirt = zeros(1,n);
h_stage = zeros(n,5);
h_interStage = zeros(1,n-1);
L = [];
for ii = 1:n

    if ii == 1
        h_aftSkirt(ii) = D(ii);                       % i-th step aft skirt length [m]
    else
        h_aftSkirt(ii) = 1/3 * D(ii) + h_dome(ii);     % i-th step aft skirt length [m]
    end

    h_interTanks(ii) = 1/4 * D(ii) + 2 * h_dome(ii);   % Inter-tank length  [m]

    if ii ~= n 
        switch AR
            case 1
                h_interStage(ii) = 5/4 * D(ii);       % Inter-stages length  [m]
            otherwise
                h_interStage(ii) = D(ii);             % Inter-stages length  [m]
        end
    else
        h_fwSkirt = 1/3 * D(ii) + h_dome(ii);         % Forward skirt length [m]
        h_interStage(ii) = h_fwSkirt;                % Forward skirt considered as the interstage with the fairing
    end
    
    h_stage(ii,:) = [h_aftSkirt(ii) h_cyl(ii,1) h_interTanks(ii) h_cyl(ii,2) h_interStage(ii)]; % Total stage length [m]
    L = [L h_stage(ii,:)];   % Lengths vector
        
end

h_PLF = 2 * D(ii);   % Payload fairing length [m]

L = [L h_PLF];      % Lengths vector

Ltot = sum(L);      % Total Launcher length

LsuD = Ltot/D(1);    % L/D of the whole launcher

h_stage_tot = sum(h_stage,2);

end