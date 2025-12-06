function [deltaV, fr0, frT, T3, guidancePara, Acs, Bcs] = EGuidance2S(rVec, V, m, mdry, mPay, thrust, Isp, guidancePara, rT, omegaT, incT, stage)
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

    rdot(cs) = dot(V,rVec)/norm(rVec);

    r(cs) = norm(rVec);
    omega(cs) = norm(cross(rVec,V))/norm(rVec)^2;

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

    lat = 90 - acosd(dot(rVec/norm(rVec), [0;0;1]));
    [~, ~, fh0, fhT] = calculateAzimuth(rVec, V, incT, lat, omegaT*rT);

    while abs((T3Old-T(3))/T(3)) > 0.001
        T3Old = T(3);
        for i = cs:2
            [r(i+1), rdot(i+1)] = verticalStateAtStaging(r(i), rdot(i), T(i), A(i), B(i), b0(i), b1(i),  c0(i), c1(i));
        end
    
        for i = cs:2
            af = a(i)/(1-T(i)/tau(i));
            fhcs0 = fh0 + (omega(i) - omega(cs))*(fhT - fh0)/(omegaT - omega(cs));
            fhcsT = fh0 + (omegaOld(i+1) - omega(cs))*(fhT - fh0)/(omegaT - omega(cs));
            omega(i+1) = horizontalStateAtStaging(A(i), B(i), T(i), omega(i), omegaOld(i+1), r(i), r(i+1), a(i), af, b0(i), b1(i), b2(i), fhcs0, fhcsT);
        end
        
        omegaOld = omega;
        
        DA = zeros([1,3]);
        DB = zeros([1,3]);
        
        for i = cs:2
            [DA(i+1), DB(i+1)] = guidanceDiscontinuities(r(i+1), rdot(i+1), omega(i+1), a(i)/(1-T(i)/tau(i)), a(i+1), Isp(i)*g0, Isp(i+1)*g0);
        end
        
        [A(cs), B(cs)] = solveEGuidanceEquations(b0, b1, c0, c1, DA, DB, T, r(cs), rT, rdot(cs));
        
        for i = cs+1:3
            A(i) = A(i-1) + DA(i) + B(i-1)*T(i-1);
            B(i) = B(i-1) + DB(i);
        end
        
        C3 = (mu/r(3)^2 - omega(3)^2*r(3))/a(3);
        if A(3) + C3 > 1
            T(3) = tau(3) - 0.01;
            break
        end

        fh30 = fh0 + (omega(3) - omega(cs))*(fhT - fh0)/(omegaT - omega(cs));
        
        T(3) = lastStageBurnTime(A(3), B(3), Isp(3), a(3), r(3), omega(3), rT, omegaT, T(3), b1(3), b2(3), fh30, fhT);
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