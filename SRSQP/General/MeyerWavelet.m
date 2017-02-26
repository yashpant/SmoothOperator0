function [phi,psi] = MeyerWavelet(t)
%% Meyer wavelet (psi) and scale (phi) as in Valenzuela and Oliviera
a1 = (2*pi)/3;
a2 = 4/3;
a3 = 4*pi/3;
a4 = pi;
a5 = (16*pi)/9;


if(t==0 || t==sqrt(a4/a5)|| t==-sqrt(a4/a5))
    %phi = 2/3 + 4/(3*pi);
    phi = (a1*cos(a1*t) + a2*cos(a3*t) - a2*a3*t*sin(a3*t))/ ...
        (a4-3*a5*t^2);
else
    phi = (sin((2*pi/3)*t) + (4*t/3)*cos(4*pi*t/3))/ ...
        (pi*t-(16*pi*t^3)/9);
end

%evalue psi_1
z=t-0.5;

c1 = 4/(3*pi);
c2 = 2*pi/3;
c3 = (1/pi);
c4 = 4*pi/3;
c5 = 16/9;

k1 = 8/(3*pi);
k2 = 8*pi/3;
k3 = 1/pi;
k4 = 4*pi/3;
k5 = 64/9;

%L'HopyBoy's rule at the points of 0/0 for psi_1
if(z==0 || z== 3/4 || z==-3/4)
    psi_1_num = -c1*c2*z*sin(c2*z)+c1*cos(c2*z)-c4*c3*cos(c4*z);
    psi_1_den = 1-3*c5*z^2;
    psi_1 = psi_1_num/psi_1_den;
    [z psi_1];
else
    psi_1_num = (4/(3*pi))*(t-0.5)*cos((2*pi/3)*(t-0.5))- ...
        (1/pi)*sin((4*pi/3)*(t-0.5));
    psi_1_den = (t-0.5)-(16/9)*(t-0.5)^3;
    psi_1 = psi_1_num/psi_1_den;
end
% for psi_2 now
if(z==0 || z==3/8 || z==-3/8)
    psi_2_num = -k1*k2*z*sin(k2*z)+k1*cos(k2*z)+k4*k3*cos(k4*z);
    psi_2_den = 1-3*k5*z^2;
    psi_2 = psi_2_num/psi_2_den;
    [z psi_2];
else
    psi_2_num = (8/(3*pi))*(t-0.5)*cos((8*pi/3)*(t-0.5))+ ...
        (1/pi)*sin((4*pi/3)*(t-0.5));
    psi_2_den = (t-0.5)-(64/9)*(t-0.5)^3;
    psi_2 = psi_2_num/psi_2_den;
end

psi = psi_1 + psi_2;
