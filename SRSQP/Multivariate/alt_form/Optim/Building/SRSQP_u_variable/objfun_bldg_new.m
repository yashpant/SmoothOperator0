function [f,g] = objfun_bldg_new(u,optParams,I1)
% maximize robustness, or min neg of robustness, bitch

xx = optParams.A_x0*optParams.x0 + optParams.B_U*u + optParams.B_D*optParams.D;
room_temp = xx(4:4:end);

[f,g_x]=getBldgRobustness(room_temp,optParams,I1,1);
f = -f;
g = -(optParams.B_U(4:4:end,:)'*g_x); %corresponding to room temp
