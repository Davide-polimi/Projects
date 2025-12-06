function Cn = Cn_computation(a,b,x,ln,d_nose,nose_type,M,alpha,h,phi,delta,A_wing,Asp_ratio_w)
d = max(a); % reference diameter assigned as the maximum semimajor axis
l = x(end);
A_ref = pi*d^2/4;

if nargin < 10
    phi = 0; % assume zero inclination of elliptical shape if not given
end

if nargin < 11
    Cn_w = 0;
end

% Cn body
Cn_ratio_SB = @(a,b) a/b*cos(phi)^2 + b/a*sin(phi)^2;
Cn_ratio_N = @(a,b) (3/2) * sqrt(a/b) * ((-b^2 / a^2) / (1 - (b^2 / a^2)^3) * log((a / b) * (1 + sqrt(1 - (b^2 / a^2))) + ...
        1 / (1 - (b^2 / a^2))));

n = max(size(a));
r = sqrt(a.*b);

dA = pi*(a(2)*b(2) - a(1)*b(1));
int_SB = (Cn_ratio_SB(a(2),b(2)) + 1) * dA/2;
int_N = Cn_ratio_SB(a(2),b(2))*r(2)*(x(2)-x(1))/2;



for i = 3:n
    dA = pi*(a(i)*b(i) - a(i-1)*b(i-1));
    int_SB = int_SB + (Cn_ratio_SB(a(i),b(i)) + Cn_ratio_SB(a(i-1),b(i-1))) * dA/2;
    if a~=b
        int_N = int_N + (Cn_ratio_N(a(i),b(i))*r(i) + Cn_ratio_N(a(i-1),b(i-1))*r(i-1))*(x(i)-x(i-1))/2;  
    end
end

if M <= 1
    eta = 0.05*(l/d)+0.52;
elseif M > 1
    eta = 1;
end

M_n = M*sin(alpha);
[~,a,~,~,nu] = atmosisa(h,extended="on");
v = M*a;
Re = v*d/nu;
Re_n = Re*sin(alpha);

%%%%%%%%%%%%%%
% ADD correlation for critical reynolds !!!!!!!!!!!!!!!!
%%%%%%%%%%%%%%


% Graphical interpolation of crossflow drag coefficient as a function of
% crossflow reynolds

Cdn_Mn = [0, 1.217;
        0.1, 1.217;
        0.2, 1.217;
        0.3, 1.223;
        0.4, 1.282;
        0.5, 1.363;
        0.6, 1.487;
        0.7, 1.537;
        0.8, 1.446;
        0.9, 1.747;
        1.0, 2.051;
        1.1, 1.819;
        1.2, 1.687;
        1.3, 1.585;
        1.4, 1.569;
        1.5, 1.545;
        1.6, 1.511;
        1.7, 1.493;
        1.8, 1.463;
        1.9, 1.439;
        2.0, 1.432;
        2.1, 1.423;
        2.2, 1.413;
        2.3, 1.399;
        2.4, 1.388;
        2.5, 1.385;
        2.6, 1.382;
        2.7, 1.375;
        2.8, 1.368;
        2.9, 1.368;
        3.0, 1.367;
        3.1, 1.362;
        3.2, 1.357;
        3.3, 1.351;
        3.4, 1.346;
        3.5, 1.341;
        3.6, 1.337;
        3.7, 1.337;
        3.8, 1.336;
        3.9, 1.330;
        4.0, 1.326;
        4.1, 1.326;
        4.2, 1.324;
        4.3, 1.317;
        4.4, 1.311;
        4.5, 1.311];

if M_n <= 4.5
    % Interpolate to find the corresponding C_dn value
    C_dn = interp1(Cdn_Mn(:,1), Cdn_Mn(:,2), M_n, 'linear');
elseif M_n > 4.5
    C_pstag = 1.8;
    C_dn = 2/3*C_pstag;
end

Cn_b = (sin(2*alpha)*cos(alpha/2)*int_SB + 2*eta*C_dn*sin(alpha)^2*int_N)/A_ref;

% Cn wing
if nargin == 13
    alpha_att = alpha + delta; % wing effective angle of attack
    if M^2 >= 1+(8/pi/Asp_ratio_w)^2
        Cn_w = (4*sin(alpha_att)*cos(alpha_att)/(M^2-1)^1/2+2*sin(alpha_att)^2)*A_wing/A_ref;
    elseif M^2 < 1+(8/pi/Asp_ratio_w)^2
        Cn_w = ((pi*Asp_ratio_w/2)*sin(alpha_att)*cos(alpha_att)+2*sin(alpha_att)^2)*A_wing/A_ref;
    end
end

Cn = Cn_b + Cn_w;

end