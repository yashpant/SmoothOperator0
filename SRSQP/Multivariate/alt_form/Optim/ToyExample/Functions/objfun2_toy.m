function [f,g] = objfun2_toy(x,optParams)
% x is the entire search vector, including inputs

%note, this guy is being minimized
%f = norm(x);
% load('OptParams.mat');



%%

g = zeros(numel(x),1);

if(optParams.robCost)
traj = reshape(x(1:optParams.dim*optParams.len),optParams.dim,optParams.len);
%f1 = alt_getRobustnessP_vector(traj,optParams.P1,optParams.Params_P1,0);

if(nargout>1)
[f1,g1] = alt_getRobustnessP_and_der_vector(traj,optParams.Params_P_unsafe);
[f2,g2]=robustness_eventually_P(traj,optParams.Params_P_term,1); %1 for der
rob = SoftMin([f1 f2]);
else
g1 = 0;g2=0;
f1 = alt_getRobustnessP_vector(traj,optParams.P_unsafe,optParams.Params_P_unsafe,0);
f2 = robustness_eventually_P(traj,optParams.Params_P_term,0);
rob = SoftMin([f1 f2]);
end

else % if no robustness in cost
f1 = 0;f2=0;rob=0;g1=0;g2=0;
end
%f2 = alt_getRobustnessP_vector(traj,optParams.P2,optParams.Params_P2,0);
%f = SoftMin([f1,f2]);
%x
%f1
f3 = (optParams.gamma)*norm(x(1:optParams.dim*optParams.len))^2; %weighted sos of states
f = -rob+f3; %min neg of robustness, robustness = not in P1
g(1:optParams.dim*optParams.len) = -g1-g2+2*x(1:optParams.dim*optParams.len)*optParams.gamma; %FIX!!!!

