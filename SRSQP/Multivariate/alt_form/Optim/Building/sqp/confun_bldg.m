function [c,ceq] = confun_bldg(u,optParams)

Nu = size(optParams.B,2);
u = reshape(u,Nu,optParams.len-1);

z = SimBldg(u,optParams);
c_state = zeros(numel(optParams.P_feas.b),optParams.len);
c_inp = zeros(numel(optParams.U_feas.b),optParams.len-1);
%State_A = kron(optParams.P_feas.A,eye(optParams.len));
%State_b = repmat(optParams.P_feas.b,optParams.len,1);
%c_state = State_A*z(:)-State_b;

%%
for i = 1:optParams.len
     c_state(:,i) = optParams.P_feas.A*z(:,i)-optParams.P_feas.b;
     if(i<=optParams.len-1)
    c_inp(:,i) = optParams.U_feas.A*u(:,i)-optParams.U_feas.b;
     end    
end

c = [c_state(:);c_inp(:)];
ceq = [];