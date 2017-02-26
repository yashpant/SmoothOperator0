function [c, ceq] = confun_linear(x,optParams)


% load('OptParams.mat');



%ceq  = x(1:2) - optParams.x0;
ceq = [];
A_big = kron(eye(optParams.len),optParams.P_feas.A);
B_big = repmat(optParams.P_feas.b,optParams.len,1);

c = A_big*x - B_big;
