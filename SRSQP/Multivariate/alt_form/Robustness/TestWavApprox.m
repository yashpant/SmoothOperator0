%%
wavparams = optParams.Params_P_unsafe;
C_00k = wavparams.C_00k;
D_ejk = wavparams.D_ejk;
k_min = wavparams.k_min;
k_max = wavparams.k_max;
j_min = wavparams.j_min;
j_max = wavparams.j_max;
E_dash = wavparams.E_dash;

%%
x = [1;0];
for i = 1:1000
clc;
tic;
alt_getWavApprox_vector_genable_spedup(x,wavparams.C_00k, ...
    wavparams.D_ejk, wavparams.k_min, wavparams.k_max, ...
    wavparams.j_min, wavparams.j_max, wavparams.E_dash);
t1(i) = toc;


tic;
alt_getWavApprox_vector_genable(x,wavparams.C_00k, ...
    wavparams.D_ejk, wavparams.k_min, wavparams.k_max, ...
    wavparams.j_min, wavparams.j_max, wavparams.E_dash);
t2(i)=toc;
end
