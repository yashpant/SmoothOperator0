function [r_P,r_P_der] = robustness_always_IfAthenB(traj,wavparamsA,wavparamsB,need_derivative)
% robustness for always_I ifAthenB, i.e. always_I (notA or B)
%trajectory Sx (n*N) (n = sys dim, N = steps)
% wavparams is a struct, {1} for A, {2} for B
%%
% if(nargin==2)
%    need_derivative = 0; 
% end
r_P_der = [];
%% hardcoded part here
pos= traj(1:3,:);
heights= traj(3,:);
wavparams1 = wavparamsA;
wavparams2 = wavparamsB;

%% the actual computation
signed_dists_inex_x = zeros(size(traj,2),1);
signed_dists_inex_z = zeros(size(traj,2),1);
smaxs = zeros(size(traj,2),1);

parfor t = 1:size(traj,2) %for all time steps
x = pos(:,t);
z = heights(:,t);
signed_dists_inex_x(t) = -alt_getWavApprox_vector_genable(x,wavparams1.C_00k, wavparams1.D_ejk,...
        wavparams1.k_min,wavparams1.k_max,wavparams1.j_min,wavparams1.j_max, ...
        wavparams1.E_dash);
signed_dists_inex_z(t) = alt_getWavApprox_vector_genable(z,wavparams2.C_00k, wavparams2.D_ejk,...
        wavparams2.k_min,wavparams2.k_max,wavparams2.j_min,wavparams2.j_max, ...
        wavparams2.E_dash);
smaxs(t) = SoftMax([signed_dists_inex_x(t) signed_dists_inex_z(t)]);    

end
Sd = smaxs; %check + - based on formula
[r_P,C] = SoftMin(Sd);


%% derivative too
