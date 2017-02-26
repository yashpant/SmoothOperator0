function [r_P,r_P_der] = alt_getRobustnessP_and_der_vector_genable(traj,wavparams,der)
% robustness for always_I notP
%trajectory Sx (n*N) (n = sys dim, N = steps)
%%

signed_dists_inex = zeros(size(traj,2),1);

parfor t = 1:size(traj,2) %for all time steps
x = traj(:,t);
signed_dists_inex(t) = alt_getWavApprox_vector_genable(x,wavparams.C_00k, wavparams.D_ejk,...
        wavparams.k_min,wavparams.k_max,wavparams.j_min,wavparams.j_max, ...
        wavparams.E_dash);
end
Sd = -signed_dists_inex; %check + - based on formula
[r_P,C] = SoftMin(Sd);


%% derivative too
%der = (nargout==2);
r_P_der = [];
if(der)
    derivative_temp = zeros(size(traj));
    for t = 1:size(traj,2) % for all steps
        dSmin_by_dut = exp(-C*Sd(t))/sum(exp(-C*Sd));
        x = traj(:,t);
        parfor i = 1:size(traj,1) %for all dims
            dSd_by_dxt_i = -alt_getWavApprox_vector_der_genable(x,i, ...
                wavparams.C_00k, wavparams.D_ejk,wavparams.k_min, ...
                wavparams.k_max, wavparams.j_min,wavparams.j_max, ...
                wavparams.E_dash);
            derivative_temp(i,t) = dSmin_by_dut*dSd_by_dxt_i;
        end
    end
    r_P_der = derivative_temp(:);
end
