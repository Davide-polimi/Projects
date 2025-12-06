function [deltaV, fr0, frT, T, A, B] = EGuidance(rVec0, VVec0, rT, VT, incT, ve, a0, T0)
    mu = 3.986004418e14;

    r0 = norm(rVec0);
    pos = rVec0/r0;
    omega0 = norm(cross(pos,VVec0))/r0;
    rdot0 = dot(pos, VVec0);

    omegaT = VT/rT;

    T = T0;
    TOld = 0;

    tau = ve/a0;

    if T > tau
        T = tau - 1;
    end
    h0 = norm(cross(rVec0, VVec0));
    hT = rT*VT;
    Dh = hT - h0;
    rm = (rT + r0)/2;
    
    lat = 90 - acosd(dot(pos, [0;0;1]));

    [~, ~, fh0, fhT] = calculateAzimuth(rVec0, VVec0, incT, lat, VT);

    deltaV = 0;
    fr0 = 0;
    frT = 0;

    i = 0;
    while abs(T-TOld)/T > 0.0001
        i = i+1;
        b0 = -ve*log(1-T/tau);
        b1 = b0*tau - ve*T;
        b2 = b1*tau - ve*T^2/2;
        c0 = b0*T - b1;
        c1 = c0*tau - ve*T^2/2;
            
        MB1 = 0 - rdot0;
        MB2 = rT - r0 - rdot0*T;

        B = (MB2 - c0*MB1/b0)/(c1 - c0*b1/b0);
        A = (MB1 - b1*B)/b0;
    
        C0 = (mu/r0^2-omega0^2*r0)/a0;

        fr0 = A + C0;
        frT = A + B*T;

        fdotr = (frT-fr0)/T;
        
        fdoth = (fhT-fh0)/T;

        ftheta = 1 - fr0^2/2 - fh0^2/2;
        fdottheta = -(fr0*fdotr + fh0*fdoth);

        fddottheta = -(fdotr^2 + fdoth^2)/2;

        %deltaV = (Dh/rm + ve*T*(fdottheta+fddottheta*tau) + fddottheta*ve*T^2/2)/(ftheta + fdottheta*tau + fddottheta*tau^2);
        deltaV = (Dh/rm - fdottheta*b1 - fddottheta*b2)/ftheta;

        TOld = T;
        T = tau*(1 - exp(-deltaV/ve));

        if T < 0 || A>1 || T >= tau || isnan(T)
            T = 100000;
            deltaV = 100000;
            fr0 = 10;
            frT = -10;
            break
        elseif i > 1000
            T = T0;
            deltaV = -ve*log(1-T/tau);
        end
    end
end