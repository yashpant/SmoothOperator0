function [phi_dash,psi_dash] = MeyerWavelet_der_sym(t)
%% Meyer wavelet (psi) and scale (phi) as in Valenzuela and Oliviera
a1 = (2*pi)/3;
a2 = 4/3;
a3 = 4*pi/3;
a4 = pi;
a5 = (16*pi)/9;
symb = 0;

if(t==0 || t==sqrt(a4/a5)|| t==-sqrt(a4/a5)) %0,3/4,-3/5
    %phi = 2/3 + 4/(3*pi); %orig form in paper for t=0
    %phi = (a1*cos(a1*t) + a2*cos(a3*t) - a2*a3*t*sin(a3*t))/ ...
    %    (a4-3*a5*t^2); % via L'Hop
    
    %for symbolic calcs of der
    if(symb) %for symbolic shite only, L'HoppyBoy here
        syms t;
        f = (sin((2*pi/3)*t) + (4*t/3)*cos(4*pi*t/3));
        g = (pi*t-(16*pi*t^3)/9);
        fdd = diff(f,2);
        gdd = diff(g,2);
        fd = diff(f,1);
        gd = diff(g,1);
        num = fdd*gd - gdd*fd;
        den = 2*gd*gd;
        phi_dash = num/den;
        clear t;
    else
        num=(32*pi*t*((4*cos((4*pi*t)/3))/3 + (2*pi*cos((2*pi*t)/3))/3 ...
            - (16*pi*t*sin((4*pi*t)/3))/9))/3 - (pi - ...
            (16*pi*t^2)/3)*((32*pi*sin((4*pi*t)/3))/9 + ...
            (4*pi^2*sin((2*pi*t)/3))/9 + ...
            (64*pi^2*t*cos((4*pi*t)/3))/27);
        den = (pi - (16*pi*t^2)/3)*(2*pi - (32*pi*t^2)/3);
        phi_dash = num/den;
    end
    
    
else %if not those NaN cases
    %phi = (sin((2*pi/3)*t) + (4*t/3)*cos(4*pi*t/3))/ ...
    %(pi*t-(16*pi*t^3)/9);
    
    if(symb)
        syms t
        phi_dash = diff(phi,1);
    else
        phi_dash = ((4*cos((4*pi*t)/3))/3 + (2*pi*cos((2*pi*t)/3))/3 ...
            - (16*pi*t*sin((4*pi*t)/3))/9)/(pi*t - (16*pi*t^3)/9) - ...
            ((pi - (16*pi*t^2)/3)*(sin((2*pi*t)/3) + ...
            (4*t*cos((4*pi*t)/3))/3))/(pi*t - (16*pi*t^3)/9)^2;
    end
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
    %psi_1_num = -c1*c2*z*sin(c2*z)+c1*cos(c2*z)-c4*c3*cos(c4*z);
    %psi_1_den = 1-3*c5*z^2;
    %psi_1 = psi_1_num/psi_1_den;
    
    if(symb)
        syms t;
        f = cos((2*pi*(t - 1/2))/3)*((1911387046407553*t)/ ...
            4503599627370496 - 1911387046407553/9007199254740992)...
            - (5734161139222659*sin((4*pi*(t - 1/2))/3)) ...
            /18014398509481984;
        
        g = (t-0.5)-(16/9)*(t-0.5)^3;
        fd = diff(f,1);
        gd = diff(g,1);
        fdd = diff(f,2);
        gdd = diff(f,2);
        num = (fdd*gd-gdd*fd); %der via L'Hop
        den = 2*gd*gd;
        psi_1_dash = num/den;
    else
        psi_1_dash = (((16*(t - 1/2)^2)/3 - 1)* ...
            ((1911387046407553*pi*sin((2*pi*(t - 1/2))/3))/ ...
            3377699720527872 - ...
            (1911387046407553*pi^2*sin((4*pi*(t - 1/2))/3))/ ...
            3377699720527872 + (4*pi^2*cos((2*pi*(t - 1/2))/3) ...
            *((1911387046407553*t)/4503599627370496 ...
            - 1911387046407553/9007199254740992))/9) - ...
            ((1911387046407553*pi*sin((2*pi*(t - 1/2))/3))/ ...
            3377699720527872 - (1911387046407553*pi^2* ...
            sin((4*pi*(t - 1/2))/3))/3377699720527872 + ...
            (4*pi^2*cos((2*pi*(t - 1/2))/3)*((1911387046407553*t)...
            /4503599627370496 - 1911387046407553/9007199254740992))...
            /9)*((1911387046407553*pi*cos((4*pi*(t - 1/2))/3))/...
            4503599627370496 - (1911387046407553* ...
            cos((2*pi*(t - 1/2))/3))/4503599627370496 + ...
            (2*pi*sin((2*pi*(t - 1/2))/3)*((1911387046407553*t)/ ...
            4503599627370496 - 1911387046407553/9007199254740992)) ...
            /3))/(((16*(t - 1/2)^2)/3 - 1)*((32*(t - 1/2)^2)/3 - 2));
    end
    
else
    if(symb)
        syms t;
    psi_1_num = (4/(3*pi))*(t-0.5)*cos((2*pi/3)*(t-0.5))- ...
    (1/pi)*sin((4*pi/3)*(t-0.5));
    psi_1_den = (t-0.5)-(16/9)*(t-0.5)^3;
    psi_1 = psi_1_num/psi_1_den;
    psi_1_dash = diff(psi_1,1);
    else
    psi_1_dash = ((1911387046407553*pi*cos((4*pi*(t - 1/2))/3))/4503599627370496 - (1911387046407553*cos((2*pi*(t - 1/2))/3))/4503599627370496 + (2*pi*sin((2*pi*(t - 1/2))/3)*((1911387046407553*t)/4503599627370496 - 1911387046407553/9007199254740992))/3)/((16*(t - 1/2)^3)/9 - t + 1/2) - (((5734161139222659*sin((4*pi*(t - 1/2))/3))/18014398509481984 - cos((2*pi*(t - 1/2))/3)*((1911387046407553*t)/4503599627370496 - 1911387046407553/9007199254740992))*((16*(t - 1/2)^2)/3 - 1))/((16*(t - 1/2)^3)/9 - t + 1/2)^2;
    end
end
% for psi_2 now
if(z==0 || z==3/8 || z==-3/8)
    %psi_2_num = -k1*k2*z*sin(k2*z)+k1*cos(k2*z)+k4*k3*cos(k4*z);
    %psi_2_den = 1-3*k5*z^2;
    %psi_2 = psi_2_num/psi_2_den;
    if(symb)
        syms t;
        f = (8/(3*pi))*(t-0.5)*cos((8*pi/3)*(t-0.5))+ ...
            (1/pi)*sin((4*pi/3)*(t-0.5));
        g = (t-0.5)-(64/9)*(t-0.5)^3;
        fd = diff(f,1);fdd = diff(f,2);
        gd = diff(g,1);gdd = diff(g,2);
        num = (fdd*gd-gdd*fd);
        dem = 2*gd*gd;
        psi_2_dash = num/den;
    else
        psi_2_dash = (((64*(t - 1/2)^2)/3 - 1)*((1911387046407553*pi^2*sin((4*pi*(t - 1/2))/3))/3377699720527872 + (1911387046407553*pi*sin((8*pi*(t - 1/2))/3))/422212465065984 + (64*pi^2*cos((8*pi*(t - 1/2))/3)*((1911387046407553*t)/2251799813685248 - 1911387046407553/4503599627370496))/9) + ((128*t)/3 - 64/3)*((1911387046407553*cos((8*pi*(t - 1/2))/3))/2251799813685248 + (1911387046407553*pi*cos((4*pi*(t - 1/2))/3))/4503599627370496 - (8*pi*sin((8*pi*(t - 1/2))/3)*((1911387046407553*t)/2251799813685248 - 1911387046407553/4503599627370496))/3))/(((16*(t - 1/2)^2)/3 - 1)*((32*(t - 1/2)^2)/3 - 2));
    end
    
else
    
    if(symb)
        syms t
        psi_2_num = (8/(3*pi))*(t-0.5)*cos((8*pi/3)*(t-0.5))+ ...
            (1/pi)*sin((4*pi/3)*(t-0.5));
        psi_2_den = (t-0.5)-(64/9)*(t-0.5)^3;
        psi_2 = psi_2_num/psi_2_den;
        psi_2_dash = diff(psi_2,1);
    else
        psi_2_dash = (((5734161139222659*sin((4*pi*(t - 1/2))/3))/18014398509481984 + cos((8*pi*(t - 1/2))/3)*((1911387046407553*t)/2251799813685248 - 1911387046407553/4503599627370496))*((64*(t - 1/2)^2)/3 - 1))/((64*(t - 1/2)^3)/9 - t + 1/2)^2 - ((1911387046407553*cos((8*pi*(t - 1/2))/3))/2251799813685248 + (1911387046407553*pi*cos((4*pi*(t - 1/2))/3))/4503599627370496 - (8*pi*sin((8*pi*(t - 1/2))/3)*((1911387046407553*t)/2251799813685248 - 1911387046407553/4503599627370496))/3)/((64*(t - 1/2)^3)/9 - t + 1/2);
    end
    
end

psi_dash = psi_1_dash + psi_2_dash;

%% arbit garbage plotting
if(0)
    ct = 0;
    for t = -5:.01:5
        ct = ct+1;
        [phi_ct(ct),psi_ct(ct)] = MeyerWavelet(t);;
        [phi_d(ct),psi_d(ct)] = MeyerWavelet_der(t);;
        clc;
        
    end
    figure;
    subplot(211)
    grid on;
    plot(-5:.01:5,phi_ct);hold all;
    plot(-5:.01:5,psi_ct); legend('phi','psi');
    subplot(212)
    grid on;
    plot(-5:.01:5,phi_d);hold all;
    plot(-5:.01:5,psi_d); legend('phid','psid');
end
