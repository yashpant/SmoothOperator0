
%%
E_dash = [0;1];
C_00k= zeros(1,1,81);
C_ejk(1,1,:) = wavparams.C(:);
k_min = wavparams.k_min;
k_max = wavparams.k_max;
j_min = wavparams.j_min;
j_max = wavparams.j_max;
D_ejk = zeros(1,3,81);
D_ejk(1,:,:) = wavparams.D;
clc;
%10*alt_getWavApprox_vector_genable(40/10,C_00k,D_ejk,k_min,k_max,j_min,j_max,E_dash)
10*getWavApprox(40/10,wavparams.C,wavparams.D,wavparams.k_min,wavparams.k_max, ...
    wavparams.j_min,wavparams.j_max)

%der


1*getWavApprox_der(40/10,wavparams.C,wavparams.D,wavparams.k_min,wavparams.k_max, ...
    wavparams.j_min,wavparams.j_max)

10*getWavApprox(40/10,wavparams.C,wavparams.D,wavparams.k_min,wavparams.k_max, ...
    wavparams.j_min,wavparams.j_max)-10*getWavApprox(39/10,wavparams.C,wavparams.D,wavparams.k_min,wavparams.k_max, ...
    wavparams.j_min,wavparams.j_max)

%%
figure;
xx = optParams.A_x0*optParams.x0 + optParams.B_U*u_0' + optParams.B_D*optParams.D;
zz = SimBldg(u_0,optParams);
plot(xx(4:4:end,:));hold all;plot(zz(4,:));

%%
figure;
xxx = optParams.A_x0*optParams.x0 + optParams.B_U*u' + optParams.B_D*optParams.D;
zzz = SimBldg(u,optParams);
plot(xxx(4:4:end,:));hold all;plot(zzz(4,:));