function [f,g] = objfun_case_u_grad(u,optParams,grad_on)

%note, this guy is being minimized
X = optParams.A_x0_nQuad*optParams.x0 + optParams.B_U_nQuad*u; %get states

g = zeros(numel(u),1);
[f1,~] = robustness_CaseStudy_u(X,optParams,grad_on);

f2 = 0;%(optParams.gamma)*norm(x(1:optParams.dim*optParams.len))^2; %weighted sos of states
f = -f1+f2; %maximize robustness, i.e. min (-robustness)
%g(1:optParams.dim*optParams.len) = -g1+2*x(1:optParams.dim*optParams.len)*optParams.gamma;

