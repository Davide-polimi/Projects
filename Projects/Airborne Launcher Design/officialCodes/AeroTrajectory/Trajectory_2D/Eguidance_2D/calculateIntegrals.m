function [b0, b1, b2, c0, c1] = calculateIntegrals(T, a, Isp)
    g0 = 9.80665;
    ve = Isp*g0;
    tau = ve./a;
    
    b0 = -ve.*log(1-T./tau);

    b1 = b0.*tau - ve.*T;

    b2 = b1.*tau - ve.*T.^2/2;

    c0 = b0.*T - b1;

    c1 = c0.*tau - ve.*T.^2/2;
end