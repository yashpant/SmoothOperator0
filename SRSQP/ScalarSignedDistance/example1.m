% example 1;

% Always Not P

ScalarSignedMain;%

%%
Ntrials = 100;
err = zeros(Ntrials,1);
clc;%close all;
for i = 1:Ntrials
    
trajec(:,i) = 2*randn(10,1);
traj = trajec(:,i);
r_exact = -getRobustnessP(traj,P,wavparams,1);
r_approx = -getRobustnessP(traj,P,wavparams,0);
err(i) = r_exact-r_approx;
err(i) = err(i)/abs(r_exact);
end
[mean((err)) std((err))]
figure;
hist(err);grid on;
