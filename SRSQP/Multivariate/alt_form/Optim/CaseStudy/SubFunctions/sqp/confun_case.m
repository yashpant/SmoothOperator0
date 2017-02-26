function [c, ceq] = confun_case(x,optParams)
%%
N = optParams.len;
Nx = optParams.dim_x;
Nu = optParams.dim_u;

x_1 = x(1:(N*Nx)+(N-1)*Nu);
x_2 = x((N*Nx)+(N-1)*Nu+1:end);
[c_1, ceq_1] = confun_SingleQuad_case(x_1,optParams.x0_1,optParams);
[c_2, ceq_2] = confun_SingleQuad_case(x_2,optParams.x0_2,optParams);

c = [c_1;c_2];
ceq = [ceq_1;ceq_2];
