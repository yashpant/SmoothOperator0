function [f,g] = objfun_linear(x,optParams)

%note, this guy is being minimized
%f = norm(x);
% load('OptParams.mat');
traj = reshape(x,optParams.dim,optParams.len);
%f1 = alt_getRobustnessP_vector(traj,optParams.P1,optParams.Params_P1,0);
%f2 = alt_getRobustnessP_vector(traj,optParams.P2,optParams.Params_P2,0);
%f = SoftMin([f1,f2]);
[f1,g] = alt_getRobustnessP_and_der_vector(traj,optParams.Params_P1);

f = f1; %min robustness, robustness = not in P1