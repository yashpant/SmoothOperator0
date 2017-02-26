
%example 1 vector, always not P
if(exist('2d_data.mat'))
    load('2d_data.mat');
else
   TestMultVarSignedDist; 
end

Ntrials = 100;
len = 10;
err = zeros(Ntrials,1);
clc;%close all;

matlabpool;
parfor i = 1:Ntrials
    i
traj = .25*randn(dim,len);

r_exact = -getRobustnessP_vector(traj,P,wavparams,1);
r_approx = -getRobustnessP_vector(traj,P,wavparams,0);
err(i) = r_exact-r_approx;
err(i) = err(i)/abs(r_exact);
end
[mean((err)) std((err))]
figure;
hist(err);grid on;