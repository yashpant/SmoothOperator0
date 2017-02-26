function f = objfun_scalar(traj,optParams)

%note, this guy is being minimized
%f = norm(x);
% load('OptParams.mat');
%traj = reshape(x,optParams.dim,optParams.len);
f1 = getRobustnessP(traj,optParams.P1,optParams.Params_P1,0);
%f2 = alt_getRobustnessP_vector(traj,optParams.P2,optParams.Params_P2,0);
%f = SoftMin([f1,f2]);
f = f1; %min robustness, robustness = not in P1
fprintf('x=%d, f=%d\n', traj(1), f);