function x_feas = bldg_getOptTraj(x,optParams,disturbances)
%%
cvx_clear;
distbs = [disturbances.d1'; disturbances.d2'; disturbances.d3'];

Nx = optParams.dim;
N = optParams.len;
Nu = size(optParams.B,2);
Nd = size(optParams.Bd,2);

cvx_begin quiet

variable z(Nx,N);
variable u(Nu,N-1);
variable dummy(10,1)

minimise sum(square(dummy-25))%0;%-sum(z(4,:))

z(:,1) == x(1:Nx);
optParams.U_feas.A*u(:,1) <= optParams.U_feas.b;
for i = 2:N
z(:,i) == optParams.A*z(:,i-1) + optParams.B*u(:,i-1) + ...
    optParams.Bd*distbs(:,i-1);
optParams.P_feas.A*z(:,i) <= optParams.P_feas.b;
optParams.U_feas.A*u(:,i-1) <= optParams.U_feas.b;
end
%z(:,N) == optParams.A*z(:,N-1) + optParams.B*u(:,N-1);

%optParams.P_final.A*z(:,N) <= optParams.P_final.b;
for i = 10:19
   dummy(i-10+1) == z(4,i); 
end

cvx_end 
clc;
x_feas.z = z;
x_feas.u = u;
x_feas.x0 = [x_feas.z(:);x_feas.u(:)];

%%
figure(1)
for i = 411:414
    subplot(i)
plot(x_feas.z(i-410,:));
grid on;
end
grid on;
figure(2)
for i = 411:414
    subplot(i)
    if(i==411)
        plot(x_feas.u);
        grid on;
    else
        j = i-411;
        stairs(distbs(j,:));
        grid on;
    end
end
grid on;