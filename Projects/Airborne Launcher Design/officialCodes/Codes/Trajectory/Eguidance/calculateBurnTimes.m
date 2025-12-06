function T = calculateBurnTimes(m, mdry, thrust, Isp)
    g0 = 9.80665;
    T = (m-mdry).*(Isp*g0)./thrust;
end