function [deltaV, fr0, frT, T3, guidancePara, Acs, Bcs] = EGuidance2S2D(r0, Vtheta0, rdot0, m, mdry, mPay, thrust, Isp, guidancePara, rT, omegaT, stage)
    A = guidancePara.A;
    B = guidancePara.B;
    T = guidancePara.T;
    omegaOld = guidancePara.omegaOld;
    
    g0 = 9.80665;
    mu = 3.986004418e14;
    
    cs = stage;
    r = zeros([1,3]);
    rdot = zeros([1,3]);
    a = zeros([1,3]);
    omega = zeros([1, 3]);

    rdot(cs) = rdot0;

    r(cs) = r0;
    omega(cs) = Vtheta0/r0;

    A(1:cs-1) = zeros([1, cs-1]);
    B(1:cs-1) = zeros([1, cs-1]);

    for i = 1:3
        a(i) = thrust(i)/(sum(m(i:3))+mPay);
    end
    
    T(1:2) = [0 0];

    T(cs:2) = calculateBurnTimes(m(cs:2), mdry(cs:2), thrust(cs:2), Isp(cs:2));

    [b0, b1, b2, c0, c1] = calculateIntegrals(T, a, Isp);

    T3Old = 0;

    tau = Isp*g0./a;

    while abs((T3Old-T(3))/T(3)) > 0.005
        T3Old = T(3);
        for i = cs:2
            [r(i+1), rdot(i+1)] = verticalStateAtStaging(r(i), rdot(i), T(i), A(i), B(i), b0(i), b1(i),  c0(i), c1(i));
        end
    
        for i = cs:2
            af = a(i)/(1-T(i)/tau(i));
            omega(i+1) = horizontalStateAtStaging2D(A(i), B(i), T(i), omega(i), omegaOld(i+1), r(i), r(i+1), a(i), af, b0(i), b1(i), b2(i));
        end
        
        omegaOld = omega;
        
        DA = zeros([1,3]);
        DB = zeros([1,3]);
    
        for i = cs:2
            [DA(i+1), DB(i+1)] = guidanceDiscontinuities(r(i+1), rdot(i+1), omega(i+1), a(i)/(1-T(i)*a(i)/(Isp(i)*g0)), a(i+1), Isp(i)*g0, Isp(i+1)*g0);
        end
        
        [A(cs), B(cs)] = solveEGuidanceEquations(b0, b1, c0, c1, DA, DB, T, r(cs), rT, rdot(cs));
        
        for i = cs+1:3
            A(i) = A(i-1) + DA(i) + B(i-1)*T(i-1);
            B(i) = B(i-1) + DB(i);
        end
        
        C3 = (mu/r(3)^2 - omega(3)^2*r(3))/a(3);
        if A(3) + C3 > 1
            %T(3) = tau(3) - 0.01;
            %break
        end

        T(3) = lastStageBurnTime2D(A(3), B(3), Isp(3), a(3), r(3), omega(3), rT, omegaT, T(3), b1(3), b2(3));
    end

    deltaV = -Isp(3)*g0*log(1-T(3)/tau(3));
    C0 = (mu/r(cs)^2-omega(cs)^2*r(cs))/a(cs);
    fr0 = A(cs) + C0;
    frT = A(3) + B(3)*T(3);
    T3 = T(3);
    Acs = A(cs);
    Bcs = B(cs);

    guidancePara.omegaOld = omega;
    guidancePara.A = A;
    guidancePara.B = B;
    guidancePara.T = T;
end