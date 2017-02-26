function [r_P] = robustness_always_SafeDist_exact(traj,d_min)
% robustness for always_I norm(p1-p2)^2>=dmin^2
%trajectory Sx (n*N) (n = sys dim, N = steps)
%%
dists_sqd = zeros(size(traj,2),1);

for t = 1:size(traj,2) %for all time steps
p1 = traj(1:3,t);
p2 = traj(4:6,t);
dists_sqd(t) = norm(p1-p2)^2;
end
Sd = dists_sqd-d_min^2; %robustness = min_I ( norm(p1-p2)^2 - d_min^2)
[r_P] = min(Sd);

