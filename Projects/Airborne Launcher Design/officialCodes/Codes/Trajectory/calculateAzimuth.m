function [az, B, fh0, fhT] = calculateAzimuth(r, V, i, lat, VT)
    pos = r/norm(r);
    Bin = asind(cosd(i)/cosd(lat));
    
    east = cross([0;0;1], pos);
    east = east/norm(east);
    north = cross(pos, east);
    north = north/norm(north);

    VTan = cross(cross(pos, V), pos);

    vy = dot(VTan, north);
    vx = dot(VTan, east);

    vxrot = VT*sind(Bin) - vx;
    vyrot = VT*cosd(Bin) - vy;

    B = atan2d(vxrot, vyrot);

    az = north*cosd(-B) + cross(pos, north)*sind(-B)+ pos*dot(pos,north)*(1-cosd(-B));
    
    Bvel = atan2d(vx, vy);
    
    fh0 = sind(Bvel - B);
    fhT = sind(Bin - B);
end