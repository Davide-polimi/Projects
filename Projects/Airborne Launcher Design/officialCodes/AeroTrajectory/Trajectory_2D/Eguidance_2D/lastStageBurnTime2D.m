function T = lastStageBurnTime2D(A, B, Isp, a, r0, omega0, rT, omegaT, T, b1, b2)
    g0 = 9.80665;
    mu = 3.986004418e14;

    ve = Isp*g0;
    tau = ve/a;

    h0 = r0^2*omega0;
    hT = rT^2*omegaT;

    Dh = hT - h0;

    rm = (r0+rT)/2;

    C0 = (mu/r0^2 - omega0^2*r0)/a;
    %CT = (mu/rT^2 - omegaT^2*rT)/a;

    fr0 = A + C0;
    TOld = 0;
    
    while abs((TOld-T)/T) > 0.001
        TOld = T;
        frT = A + B*T;% + CT;
    
        fdr = (frT - fr0)/T;
    
        ft = 1 - fr0^2/2;
        fdt = -(fr0*fdr);
        fddt = -(fdr^2)/2;
    
        DV = (Dh/rm - fdt*b1 - fddt*b2)/ft;
    
        T = tau*(1 - exp(-DV/ve));
    end
end