function [r_P] = robustness_eventually_P_exact(traj,Set)
% robustness for eventually_I inP
%trajectory Sx (n*N) (n = sys dim, N = steps)
%%

signed_dists_inex = zeros(size(traj,2),1);

for t = 1:size(traj,2) %for all time steps
x = traj(:,t);
signed_dists_inex(t) = SignedDist(x,Set.A,Set.b);
end
Sd = signed_dists_inex; %check + - based on formula
r_P = max(Sd);



