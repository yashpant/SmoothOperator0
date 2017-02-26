function [phi_dash] = MeyerWavelet_der_phi(t)
%% Meyer wavelet (psi) and scale (phi) as in Valenzuela and Oliviera


if(t==0 || t==3/4 || t==-3/4) %0,3/4,-3/5
    %phi = 2/3 + 4/(3*pi); %orig form in paper for t=0
    %phi = (a1*cos(a1*t) + a2*cos(a3*t) - a2*a3*t*sin(a3*t))/ ...
    %    (a4-3*a5*t^2); % via L'Hop
    
    
        num=(32*pi*t*((4*cos((4*pi*t)/3))/3 + (2*pi*cos((2*pi*t)/3))/3 ...
            - (16*pi*t*sin((4*pi*t)/3))/9))/3 - (pi - ...
            (16*pi*t^2)/3)*((32*pi*sin((4*pi*t)/3))/9 + ...
            (4*pi^2*sin((2*pi*t)/3))/9 + ...
            (64*pi^2*t*cos((4*pi*t)/3))/27);
        den = (pi - (16*pi*t^2)/3)*(2*pi - (32*pi*t^2)/3);
        phi_dash = num/den;
    
    
    
else %if not those NaN cases
    %phi = (sin((2*pi/3)*t) + (4*t/3)*cos(4*pi*t/3))/ ...
    %(pi*t-(16*pi*t^3)/9);
    
    
        phi_dash = ((4*cos((4*pi*t)/3))/3 + (2*pi*cos((2*pi*t)/3))/3 ...
            - (16*pi*t*sin((4*pi*t)/3))/9)/(pi*t - (16*pi*t^3)/9) - ...
            ((pi - (16*pi*t^2)/3)*(sin((2*pi*t)/3) + ...
            (4*t*cos((4*pi*t)/3))/3))/(pi*t - (16*pi*t^3)/9)^2;
    
end

