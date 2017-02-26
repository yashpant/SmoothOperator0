function [c,ceq] = confun_toy_falsification(x0,optParams)
ceq = [];

c = optParams.P_init.A*x0-optParams.P_init.b;
