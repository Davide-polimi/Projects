function [Cl,Cd] = aero_conversion(Cn,Ca,alpha)
    % To obtain lift and drag coefficients from normal and axial coefficients

    Cl = Cn*cos(alpha)-Ca*sin(alpha);
    Cd = Cn*sin(alpha)+Ca*cos(alpha);
end