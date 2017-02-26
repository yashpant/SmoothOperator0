function [f,g] = objfun_case(x,optParams)

%note, this guy is being minimized
%f = norm(x);
% load('OptParams.mat');


g = zeros(numel(x),1);

%if(optParams.robCost)

%f1 = alt_getRobustnessP_vector(traj,optParams.P1,optParams.Params_P1,0);
[f1,~] = robustness_CaseStudy(x,optParams);
%else
%f1 = 0;
%end
f2 = 0;%(optParams.gamma)*norm(x(1:optParams.dim*optParams.len))^2; %weighted sos of states
f = -f1+f2; %maximize robustness, i.e. min (-robustness)
%g(1:optParams.dim*optParams.len) = -g1+2*x(1:optParams.dim*optParams.len)*optParams.gamma;

