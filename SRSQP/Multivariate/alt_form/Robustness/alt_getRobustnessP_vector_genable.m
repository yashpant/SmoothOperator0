function [r_P,Sd] = alt_getRobustnessP_vector_genable(traj,P,wavparams,exact)
%trajectory Sx (n*N) (n = sys dim, N = steps)
% exact: 0 = approximation

Sd = zeros(size(traj,2),1);
for i = 1:size(traj,2)
    x = traj(:,i);
    Sd(i) = -alt_getWavApprox_vector_genable(x,wavparams.C_00k, wavparams.D_ejk,...
        wavparams.k_min,wavparams.k_max,wavparams.j_min,wavparams.j_max, ...
        wavparams.E_dash);  
end
r_P = SoftMin(Sd);
