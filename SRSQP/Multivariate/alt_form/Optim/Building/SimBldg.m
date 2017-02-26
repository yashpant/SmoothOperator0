function z = SimBldg(u,optParams)

%%
distbs = [optParams.disturbances.d1'; ...
            optParams.disturbances.d2'; ...
            optParams.disturbances.d3'];
%%
Nx = optParams.dim;
N = optParams.len;

%%
z = zeros(Nx,N);

z(:,1) = optParams.x0;
for i = 2:N
z(:,i) = optParams.A*z(:,i-1) + optParams.B*u(:,i-1) + ...
    optParams.Bd*distbs(:,i-1);
end
