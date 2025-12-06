function [gridL, gridD] = generateAeroTable(rocketDesign,engine_mode,phi,wing_design, launchHeight, launchVelocity)
    V = [launchVelocity-50:10:500, 530:30:5000];
    h = (launchHeight-1000):500:84000;
    AOA = [0:0.1:30];

    [x, y, z] = ndgrid(V, h, AOA);

    tableL = zeros([numel(V), numel(h), numel(AOA)]);
    tableD = tableL;

    parfor i = 1:numel(V)
        tempTableL = zeros([numel(h), numel(AOA)]);
        tempTableD = zeros([numel(h), numel(AOA)]);
        for j = 1:size(h,2)
            for k = 1:size(AOA,2)
                [L,D] = Aerodynamics_computation(rocketDesign,V(i),AOA(k)*pi/180,h(j),engine_mode,phi,wing_design);
                tempTableL(j,k) = L/V(i)^2;
                tempTableD(j,k) = D/V(i)^2;
            end
        end
        tableL(i, :, :) = tempTableL;
        tableD(i, :, :) = tempTableD;
    end

    gridL = griddedInterpolant(x, y, z, tableL, 'linear');
    gridD = griddedInterpolant(x, y, z, tableD, 'linear');
end