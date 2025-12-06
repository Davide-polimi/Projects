function [dVdt, L, D] = calcAccel(r, V, m, thrust, pitch, incT, VT, aero)
    wE = [0; 0; 7.2921159e-5]; % Angular velocity of earth
    mu = 3.986004418e14; % Gravitational constant of earth
    rE = 6371008.7714; % Radius of earth

    lat = 90 - acosd(dot(r/norm(r), [0;0;1]));
    [az, ~] = calculateAzimuth(r, V, incT, lat, VT);
    auxVec = cross(az, r);
    auxVec = auxVec/norm(auxVec);
    aDir = az*cosd(pitch) + cross(auxVec, az)*sind(pitch)+ auxVec*dot(auxVec,az)*(1-cosd(pitch));

    L = 0;
    D = 0;
    airV = [1;0;0];
    Ldir = [1;0;0];

    if norm(r)-rE < 84000
        airV = V - cross(wE, r);
        auxAirV = cross(airV, r);
        Ldir = cross(auxAirV, airV);
        Ldir = Ldir/norm(Ldir);
        
        AOA = real(acosd(dot(aDir,airV)/norm(airV)));
        if dot(airV,r)/norm(airV) > dot(aDir,r)
            Ldir = -Ldir;
        end
        L = aero.L(norm(airV), norm(r)-rE, AOA)*norm(airV)^2;
        D = aero.D(norm(airV), norm(r)-rE, AOA)*norm(airV)^2;
    end

    dVdt = (aDir*thrust - D*airV/norm(airV) + L*Ldir)/m - mu/norm(r)^3*r;
end