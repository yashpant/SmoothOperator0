function [f,g] = objfun_bldg_new_rh(u,optParams,I1,i,x_upto)
% maximize robustness, or min neg of robustness, bitch

%xx = optParams.A_x0*optParams.x0 + optParams.B_U*u + optParams.B_D*optParams.D;

xx = x_upto;

for t = i+1:optParams.len
        
   xx(:,t) = optParams.A*xx(:,t-1) + optParams.B*u(t-1,1) ;
    
end
room_temp = xx(4,:);
[f,g_x]=getBldgRobustness(room_temp,optParams,I1,1);
f = -f;
g = -(optParams.B_U(4:4:end,:)'*g_x); %corresponding to room temp
