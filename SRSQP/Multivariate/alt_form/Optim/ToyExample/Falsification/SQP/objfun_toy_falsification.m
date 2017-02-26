function f = objfun_toy_falsification(x0,optParams,exact)

%%
len = optParams.len;
x_star = optParams.x_star;
K = optParams.K;
A = optParams.A;
B = optParams.B;
x = zeros(optParams.dim,len);
x(:,1) = x0;
for t = 2:len
    x(:,t) = (A-B*K)*x(:,t-1) + B*K*x_star;%(A-B*K)*(x(:,t-1) - x_star) + x_star;      
end
traj = x;

%f1 = alt_getRobustnessP_vector(traj,optParams.P1,optParams.Params_P1,0);
if(~exact)
f = robustness_toy_mex(traj,optParams.Params_P1);
else
f = alt_getRobustnessP_vector(x,optParams.P1,optParams.Params_P1,1);
end
