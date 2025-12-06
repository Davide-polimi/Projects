clc
clear
close all

% INPUTS
Mach = 0:0.2:8;
alpha = deg2rad(0:1:20);
h = 10000;

x = [0 4.5 21];
a = [0 1.7 1.7];
b = a;
[~,a_sound,~,~,~] = atmosisa(h);
v = Mach*a_sound;

ln = x(2);
d_nose = a(2);
nose_type = 'C';
rocket_design.a = a;
rocket_design.b = b;
rocket_design.x = x;
rocket_design.ln = ln;
rocket_design.d_nose = d_nose;
rocket_design.nose_type = nose_type;
wing_design = 0;
phi = 0;
engine_mode = 1;


for i = 1:length(Mach)
    for j = 1:length(alpha)
        [Ca(i,j),Ca_w(i,j),Ca_f(i,j),Ca_b(i,j)] = Ca_computation(a,b,x,ln,d_nose,nose_type,Mach(i),alpha(j),h,0);
        Cn(i,j) = Cn_computation(a,b,x,ln,d_nose,nose_type,Mach(i),alpha(j),h);
        [Cl(i,j),Cd(i,j)] = aero_conversion(Cn(i,j),Ca(i,j),alpha(j));    
        [L(i,j), D(i,j)] = Aerodynamics_computation(rocket_design, v(i), alpha(j), h, engine_mode, phi, wing_design);
    end
end

figure(1)
plot(x,a/2,'k',x,-a/2,'k')
axis equal
grid on
title('Launcher shape')
xlabel('Axis coordinate [m]')
ylabel('Cross dimension [m]')

figure()
plot(Mach,Ca)
xlabel('Mach number')
ylabel('CA')
grid on
title('Axial coefficient')

figure()
plot(Mach,Ca_w)
xlabel('Mach number')
ylabel('CA_w')
grid on
title('Wave drag coefficient')

figure()
plot(Mach,Ca_b)
xlabel('Mach number')
ylabel('CA_b')
grid on
title('Base drag coefficient')

figure()
plot(Mach,Ca_f)
xlabel('Mach number')
ylabel('CA_f')
grid on
title('Skin friction coefficient')

figure()
plot(rad2deg(alpha),Cn)
xlabel('Alpha')
ylabel('Cn')
grid on
title('Normal coefficient')

figure()
plot(Mach,Cd)
xlabel('Mach')
ylabel('Cd')
grid on
title('Cd(mach;alpha)')

figure()
plot(rad2deg(alpha),Cl)
xlabel('Alpha')
ylabel('Cl')
grid on
title('Cl(alpha;Mach)')


