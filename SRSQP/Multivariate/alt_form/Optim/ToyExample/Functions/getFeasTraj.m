function x_feas = getFeasTraj(x,optParams)


Nx = optParams.dim;
N = optParams.len;
Nu = size(optParams.B,2);

cvx_begin quiet

variable z(Nx,N);
variable u(Nu,N-1);

minimise 0

z(:,1) == x(1:Nx);
optParams.U_feas.A*u(:,1) <= optParams.U_feas.b;
for i = 2:N-1
z(:,i) == optParams.A*z(:,i-1) + optParams.B*u(:,i-1);
optParams.P_feas.A*z(:,i) <= optParams.P_feas.b;
optParams.U_feas.A*u(:,i) <= optParams.U_feas.b;
end
z(:,N) == optParams.A*z(:,N-1) + optParams.B*u(:,N-1);

optParams.P_final.A*z(:,N) <= optParams.P_final.b;

cvx_end 

x_feas.z = z;
x_feas.u = u;
x_feas.x0 = [x_feas.z(:);x_feas.u(:)];
