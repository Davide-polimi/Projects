function T = lastStageBurnTime(A, B, Isp, a, r0, omega0, rT, omegaT, T, b1, b2, fh0, fhT)
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
    
    while abs((TOld-T)/T) > 0.0001
        TOld = T;
        frT = A + B*T;% + CT;
    
        fdr = (frT - fr0)/T;
        fdh = (fhT - fh0)/T;
    
        ft = 1 - fr0^2/2 - fh0^2/2;
        fdt = -(fr0*fdr+fh0*fdh);
        fddt = -(fdr^2+fdh^2)/2;

        %DV = (Dh/rm + ve*T*(fdt + fddt*tau) + fddt*ve*T^2/2)/(ft + fdt*tau + fddt*tau^2);
        DV = (Dh/rm - fdt*b1 - fddt*b2)/ft;
    
        T = tau*(1 - exp(-DV/ve));

        b1 = DV.*tau - ve.*T;
        b2 = b1.*tau - ve.*T.^2/2;
    end
end