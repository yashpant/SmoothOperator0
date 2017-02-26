x0_reshape = reshape(x_0(1:optParams.len*optParams.dim),optParams.dim,optParams.len);
x_reshape = reshape(x(1:optParams.len*optParams.dim),optParams.dim,optParams.len);

figure;
plot(optParams.P_feas,'Color','gray');
grid on;
hold on;
plot(optParams.P_final,'Color','white');
hold on;
plot(optParams.P1,'Color','red');
hold on;

for i = 1:optParams.len
    plot(x0_reshape(1,i),x0_reshape(2,i),'o');
    hold on;
    plot(x_reshape(1,i),x_reshape(2,i),'g*');
    pause(0.01);
    if(i==1)
       hold on;
       plot(x0_reshape(1,i),x0_reshape(2,i),'m+');
       legend('FeasSet','TermSet','UnsafeSet','init traj','computed','x_0'); 
    end
    
end

traj = reshape(x(1:optParams.dim*optParams.len),optParams.dim,optParams.len);
exact_robustness = alt_getRobustnessP_toy(traj,optParams.P1,optParams.Params_P1,1)
smooth_robustness = alt_getRobustnessP_toy(traj,optParams.P1,optParams.Params_P1,0)