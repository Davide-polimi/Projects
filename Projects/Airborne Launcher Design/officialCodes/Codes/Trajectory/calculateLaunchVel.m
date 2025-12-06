function vel = calculateLaunchVel(launchVelocity, launchPitch, launchAzimuth, r)
    launchAzimuth = -launchAzimuth;
    pos = r/norm(r);
    vec = -cross([0;1;0], pos);
    vec = vec*cosd(launchAzimuth) + cross(pos, vec)*sind(launchAzimuth)+ pos*dot(pos,vec)*(1-cosd(launchAzimuth));
    aux = cross(vec,pos);
    vec = vec*cosd(launchPitch) + cross(aux, vec)*sind(launchPitch)+ aux*dot(aux,vec)*(1-cosd(launchPitch));
    vel = vec*launchVelocity/norm(vec);
end