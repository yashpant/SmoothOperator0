function f = objfun_bldg(u,optParams,I1)
% maximize robustness, or min neg of robustness, bitch
z = SimBldg(u,optParams);
room_temp = z(4,:);

if(optParams.exact)
   f =  -getBldgRobustness_exact(room_temp,optParams,I1);
else
   f =  -getBldgRobustness(room_temp,optParams,I1);
end