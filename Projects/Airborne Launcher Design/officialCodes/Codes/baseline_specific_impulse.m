clear
close all
clc

Isp_1stage = [330 290.2 290 279 280 335]; % s 
T_1stage = [327 726 556 4323 1371 3015 ]; % kN 
tburn_1stage = [180 68.6 169 135.7 109.9 155]; % s 

c_tburn = polyfit(Isp_1stage,tburn_1stage, 1);

Isp_1stage_vec = linspace(279,432,100);
tburn_vec = polyval(c_tburn,Isp_1stage_vec,length(Isp_1stage_vec));

figure
plot(Isp_1stage_vec,tburn_vec, LineWidth=2)
hold on
for i = 1: length(Isp_1stage)

    plot(Isp_1stage(i),tburn_1stage(i), 'r.', LineWidth=2, MarkerSize=20)

end
grid on
xlabel('Specific Impulse first stage[s]')
ylabel('Burning time first stage[s]')

Isp_2stage = [330 289.4 327 293.5 458 343 287.5 348 421];
t_burn_2stage = [360 69.4 378 92.9	900	373	77.1 397 360];
c_tburn2 = polyfit(Isp_2stage,t_burn_2stage, 1);
Isp_2stage_vec = linspace(287.5,458,100);
tburn2_vec = polyval(c_tburn2,Isp_2stage_vec, length(Isp_2stage_vec));

figure
plot(Isp_2stage_vec,tburn2_vec, LineWidth=2)
hold on
for i = 1: length(Isp_2stage)

    plot(Isp_2stage(i),t_burn_2stage(i), 'r.', LineWidth=2, MarkerSize=20)

end
grid on
xlabel('Specific Impulse second stage[s]')
ylabel('Burning time second stage[s]')
