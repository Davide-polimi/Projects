function cdw = wvdrgogv(diameter,length,Mach)

% Returns the supersonic wave drag for a tangent ogive
% by empirical method of Miles
%
% diameter = ogive base diameter
% length = lengthof ogive
% Mach = Mach number

% The values resulting from these calculations are valid in the 1.5 to 3.5
% Mach range and for semi-vertex angles between 10 and 25 degrees

if Mach > 3.5
    Mach = 3.5; % modifying Mach number into the range of model validity
end

sigma = 2*(180/pi)*atan(diameter/2/length);
P = (0.083+0.096/Mach^2)*(sigma/10)^1.69;
lod2 = (length/diameter)^2;
cdw = P*(1-(196*lod2-16)/(14*(Mach+18)*lod2));