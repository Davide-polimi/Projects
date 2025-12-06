function omega = horizontalStateAtStaging(A, B, T, omega0, omegaT, r0, rT, a0, aT, b0, b1, b2, fh0, fhT)
    mu = 3.986004418e14;
    
    omegaOld = -10000000;
    omega = omegaT;

    while abs((omegaOld - omega)/omegaOld) > 0.01
        fr0 = A + (mu/r0^2-omega0^2*r0)/a0;
        frT = A + B*T + (mu/rT^2-omega^2*rT)/aT;
    
        fdr = (frT-fr0)/T;
        fdh = (fhT-fh0)/T;
    
        ft = 1 - fr0^2/2 - fh0^2/2;
        fdt = -(fr0*fdr + fh0*fdh);
        fddt = -(fdr^2+fdh^2)/2;
    
        hT = (r0 + rT)/2*(ft*b0 + fdt*b1 + fddt*b2) + omega0*r0^2;
    
        omegaOld = omega;
        omega = hT/rT^2;
    end
end