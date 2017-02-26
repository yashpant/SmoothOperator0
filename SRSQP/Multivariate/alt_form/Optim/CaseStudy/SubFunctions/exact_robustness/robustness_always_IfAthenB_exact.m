function [r_P] = robustness_always_IfAthenB_exact(traj,SetA,SetB)
% robustness for always_I ifAthenB, i.e. always_I (notA or B)
%trajectory Sx (n*N) (n = sys dim, N = steps)
% wavparams is a struct, {1} for A, {2} for B
%%

%% hardcoded part here
pos= traj(1:3,:);
heights= traj(3,:);


%% the actual computation
signed_dists_inex_x = zeros(size(traj,2),1);
signed_dists_inex_z = zeros(size(traj,2),1);
smaxs = zeros(size(traj,2),1);

for t = 1:size(traj,2) %for all time steps
x = pos(:,t);
z = heights(:,t);
signed_dists_inex_x(t) = -SignedDist(x,SetA.A,SetA.b);
signed_dists_inex_z(t) = SignedDist(z,SetB.A,SetB.b);
smaxs(t) = max([signed_dists_inex_x(t) signed_dists_inex_z(t)]);    

end
Sd = smaxs; %check + - based on formula
[r_P] = min(Sd);


