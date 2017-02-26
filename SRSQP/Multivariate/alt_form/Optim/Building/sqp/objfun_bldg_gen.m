function f = objfun_bldg_gen(u,optParams,I1)
% maximize robustness, or min neg of robustness, bitch
z = SimBldg(u,optParams);
room_temp_1 = z(4,:);

xx = optParams.A_x0*optParams.x0 + optParams.B_U*u' + optParams.B_D*optParams.D;
room_temp = xx(4:4:end);
max(abs(room_temp-room_temp_1'))
f =  -getBldgRobustness(room_temp,optParams,I1,0);
