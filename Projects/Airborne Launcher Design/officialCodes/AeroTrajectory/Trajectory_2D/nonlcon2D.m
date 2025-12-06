function [cneq, ceq] = nonlcon2D(x, Dt, Tf, r0, V0, m0, mDry, mPayload, rT, VT, Isp, thrust, launchLatitude, rocketDesign, nStages)
    g0 = 9.80665; % g0
    rE = 6371008.7714; % Radius of earth
    wE = 7.2921159e-5; % Angular velocity of earth
    mu = 3.986004418e14; % Gravitational constant of earth

    cs = 4 - nStages;

    ve = Isp(cs) * g0;

    ySol = [r0; V0'; sum(m0) + mPayload]; % [r, Vtheta, rdot, m]
    qMax = 0;

    for i = 1:floor(Tf/Dt)
        pitch0 = x(i);
        pitchT = x(i+1);
        pitch = @(t) pitch0 - (pitch0 - pitchT)*t/Dt;

        airV = [ySol(2)+wE*ySol(1)*cosd(launchLatitude),ySol(3)];
        [~,~,~,rho] = atmosisa(ySol(1)-rE,'extended', 'on');
        q = 1/2*rho*norm(airV)^2;
        gAirV = atan2d(airV(2),airV(1));
        
        AOA = (pitch0+pitchT)/2 - gAirV;
        
        [L,D] = Aerodynamics_computation(rocketDesign,norm(airV),abs(AOA)*pi/180,ySol(1)-rE,1,0,0);

        ar = @(t,y) (sind(pitch(t))*thrust(cs) - D*sind(gAirV) + L*cosd(gAirV))/y(4) - mu/y(1)^2 + y(2)^2/y(1);
        at = @(t,y) (cosd(pitch(t))*thrust(cs) - D*cosd(gAirV) - L*sind(gAirV))/y(4) - 2*y(3)*(y(2)/y(1));

        yd = @(t,y) [y(3); ar(t,y); at(t,y); -thrust(cs)/ve];
        
        [~, ySol] = ode45(yd, [0 Dt], ySol);
        ySol = ySol(end,:)';
        
        if q > qMax
            qMax = q;
        elseif q < 0.05*qMax
            break
        end
    end

    r0 = ySol(1);
    V0 = ySol(2:3);

    m = m0;
    m(cs) = ySol(4) - sum(m0(cs+1:3)) - mPayload;

    guidancePara.omegaOld = [0.005 0.006 0.007];
    guidancePara.A = [0.5 0.5 0.5];
    guidancePara.B = [0 0 0];
    guidancePara.T = [0 0 (m0(3)-mDry(3))*Isp(3)*g0/thrust(3)];

    [~, fr0, frT] = EGuidance2S2D(r0, V0(1), V0(2), m, mDry, mPayload, thrust, Isp, guidancePara, rT, VT/rT, cs);

    cneq = [fr0 - 1, -frT - 1];

    ceq = [];
end