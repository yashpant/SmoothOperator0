function [r_P,r_P_der] = robustness_eventually_P_genable(traj,wavparams,need_derivative)
% robustness for eventually_I inP
%trajectory Sx (n*N) (n = sys dim, N = steps)

%%
signed_dists_inex = zeros(size(traj,2),1);

parfor t = 1:size(traj,2) %for all time steps
x = traj(:,t);
signed_dists_inex(t) = alt_getWavApprox_vector_genable(x,wavparams.C_00k, wavparams.D_ejk,...
        wavparams.k_min,wavparams.k_max,wavparams.j_min,wavparams.j_max, ...
        wavparams.E_dash);
end
Sd = signed_dists_inex; %check + - based on formula
[r_P,C] = SoftMax(Sd);


%% derivative too
% der = need_derivative; %to be replaced by nargout
% derivative_temp = zeros(size(traj));

der = need_derivative;r_P_der = [];
if(der)
    %derivative_temp = zeros(size(traj));
    for t = 1:size(traj,2) % for all steps
        dSmax_by_dut = exp(C*Sd(t))/sum(exp(C*Sd));
        x = traj(:,t);
        parfor i = 1:size(traj,1) %for all dims
            dSd_by_dxt_i = alt_getWavApprox_vector_der_genable(x,i, ...
                wavparams.C_00k, wavparams.D_ejk,wavparams.k_min, ...
                wavparams.k_max, wavparams.j_min,wavparams.j_max, ...
                wavparams.E_dash);
            derivative_temp(i,t) = dSmax_by_dut*dSd_by_dxt_i;
        end
    end
    r_P_der = derivative_temp(:);
end
%end
