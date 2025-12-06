function [DA, DB] = guidanceDiscontinuities(rT, rdT, omegaT, a0T, a1T, ve0, ve1)
    mu = 3.986004418e14;

    DA = (mu/rT^2 - omegaT^2*rT)*(1/a0T - 1/a1T);
    DB = (mu/rT^2 - omegaT^2*rT)*(1/ve0 - 1/ve1) + (3*omegaT^2 - 2*mu/rT^3)*rdT*(1/a0T - 1/a1T);
end