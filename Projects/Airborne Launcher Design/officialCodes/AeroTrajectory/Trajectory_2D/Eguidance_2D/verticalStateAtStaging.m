function [r, rdot] = verticalStateAtStaging(r0, rdot0, T, A, B, b0, b1, c0, c1)
    rdot = rdot0 + b0*A + b1*B;
    r = r0 + rdot0*T + c0*A + c1*B;
end