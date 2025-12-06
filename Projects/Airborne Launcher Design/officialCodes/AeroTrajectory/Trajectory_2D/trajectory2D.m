function [qMaxLoad, deltaV] = trajectory2D(Dt, structure, propulsion, launcher, stability, rocketDesign, target)
    
    % Phisics constants
    g0 = 9.80665; % g0
    rE = 6371008.7714; % Radius of earth
    wE = 7.2921159e-5; % Angular velocity of earth
    mu = 3.986004418e14;
    
    % Design constants
    mDry = structure.mDry; % [kg] Dry mass
    m0 = structure.m0; % [kg] Initial mass
    mPayload = structure.mPayload; % [kg] Payload mass
    Isp = propulsion.Isp; % [s] Isp
    thrust = propulsion.thrust; % [N] Thrust
    launchVelocity = launcher.launchVelocity; % [m/s]
    launchPitch = launcher.launchPitch; % [degrees]
    launchHeight = launcher.launchHeight;
    launchLatitude = launcher.launchLatitude; % [degrees]
    pitchRate = stability.pitchRate; % [degrees/second]
    pitchAccel = stability.pitchAccel; % [degrees/second^2]
    
    % Target Orbit
    rT = target.r;
    VT = target.V;
    
    % Initial Values
    r0 = rE + launchHeight;
    V0 = [launchVelocity*cosd(launchPitch) - wE*r0*cosd(launchLatitude), launchVelocity*sind(launchPitch)];
    
    % Simulation constants
    Tf = (m0(2)-mDry(2))*Isp(2)*g0/thrust(2);
    nPoints = floor(Tf/Dt) + 1;
    
    x0 = [launchPitch:pitchRate*Dt:40, ones([1, nPoints])*40]'; % [Pitch]
    x0 = x0(1:nPoints);
    
    % Setup linear constraints
    Aineq = zeros([(nPoints-1)*4, nPoints]);
    bineq = zeros([(nPoints-1)*4, 1]);
    
    for i = 0:nPoints-2
        % Max pitch rate
        Aineq(4*i+1,:) = [zeros([1, i]), 1, -1, zeros([1, nPoints-2-i])];
        Aineq(4*i+2,:) = -Aineq(4*i+1,:);
        bineq(4*i+1) = pitchRate*Dt;
        bineq(4*i+2) = pitchRate*Dt;
        if i < nPoints-2 % Max pitch acceleration
            Aineq(4*i+3,:) = [zeros([1, i]), 1, -2, 1, zeros([1, nPoints-3-i])];
            Aineq(4*i+4,:) = -Aineq(4*i+3,:);
            bineq(4*i+3) = pitchAccel*Dt^2;
            bineq(4*i+4) = pitchAccel*Dt^2;
        end
    end
    % Start with 0 pitch rate
    Aineq(end-1,:) = [-1 1 zeros([1, nPoints-2])];
    Aineq(end,:) = -Aineq(end-1,:);
    bineq(end-1) = pitchAccel*Dt^2;
    bineq(end) = pitchAccel*Dt^2;
    
    % Initial and final pitch
    Aeq = [1, zeros([1, nPoints-1])];
    beq = [launchPitch];
    
    %%
    problem.x0 = x0;
    problem.Aineq = Aineq;
    problem.bineq = bineq;
    problem.Aeq = Aeq;
    problem.beq = beq;
    problem.lb = -ones([nPoints, 1])*90;
    problem.ub = ones([nPoints, 1])*90;
    problem.solver = 'fmincon';
    problem.options = optimoptions('fmincon','Algorithm', 'sqp','Display','iter',...
        'UseParallel',true,'MaxFunctionEvaluations', 3000);
    problem.objective = @(x) objective2D(x, Dt, Tf, r0, V0, m0, mDry, mPayload, rT, VT, Isp, thrust, launchLatitude, rocketDesign, 2);
    problem.nonlcon = [];%@(x) nonlcon2D(x, Dt, Tf, r0, V0, m0, mDry, mPayload, rT, VT, Isp, thrust, launchLatitude, rocketDesign, 2);
    
    x = fmincon(problem);
    
    %%
    cs = 2;
    
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
    
        yd = @(t,y) [y(3); at(t,y); ar(t,y); -thrust(cs)/ve];
        
        [~, ySol] = ode45(yd, [0 Dt], ySol);
        ySol = ySol(end,:)';
        
        if q > qMax
            qMax = q;
            qMaxLoad.q = qMax;
            qMaxLoad.L = L;
            qMaxLoad.D = D;
            qMaxLoad.a = thrust(cs)/ySol(4);
            qMaxLoad.AOA = AOA;
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
    [deltaV2] = EGuidance2S2D(r0, V0(1), V0(2), m, mDry, mPayload, thrust, Isp, guidancePara, rT, VT/rT, cs);
    deltaV = Isp(1)*g0*log((sum(m0)+mPayload)/(mDry(1)+m0(2)+m0(3)+mPayload)) + Isp(2)*g0*log((m0(2)+m0(3)+mPayload)/(mDry(2)+m0(3)+mPayload)) + deltaV2;
end