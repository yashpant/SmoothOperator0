function [r_notP] = robustness_toy_genable(traj,wavparams)
% robustness for always_I notP
%trajectory Sx (n*N) (n = sys dim, N = steps)
%%

signed_dists_inex = zeros(size(traj,2),1);

for t = 1:size(traj,2) %for all time steps
x = traj(:,t);
signed_dists_inex(t) = alt_getWavApprox_vector_genable(x,wavparams.C_00k, wavparams.D_ejk,...
        wavparams.k_min,wavparams.k_max,wavparams.j_min,wavparams.j_max, ...
        wavparams.E_dash);
end
Sd = -signed_dists_inex; %check + - based on formula
[r_notP] = SoftMin(Sd);