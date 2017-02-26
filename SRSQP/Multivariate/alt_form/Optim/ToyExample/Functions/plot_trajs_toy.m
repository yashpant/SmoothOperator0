
Nx = optParams.dim;
N = optParams.len;
Nu = size(optParams.B,2);

traj_0 = reshape(x_0(1:Nx*N),Nx,N);
traj = reshape(x(1:Nx*N),Nx,N);

figure;
plot(optParams.P_feas,'color','green');
hold on;
plot(optParams.P1);
hold on;grid on;
plot(optParams.P_final,'color','gray');
for i = 1:N
    hold on;
    plot(traj(1,i),traj(2,i),'k*');
    hold on
    plot(traj_0(1,i),traj(2,i),'bo');
end