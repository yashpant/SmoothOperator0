
%example 1 vector, always not P
if(exist('newdata_R3.mat'))
    load('newdata_R3.mat');i
	'loading data'
	
else
   %alt_TestMultVarSignedDist; 
	'file no exist, no'
	keyboard
end

Ntrials = 50;
len = 250;
err = zeros(Ntrials,1);
clc;%close all;

x_lb = -2.5;
x_ub = 2.5;
u_lb = -0.32;
u_ub = 0.32;
Nx = 3;
Nu = 3;
A = eye(3);
B = eye(3);
traj = zeros(Nx,len);
Traj_all = zeros(Nx,len,Ntrials);
U_all = zeros(Nx,len-1,Ntrials);
%matlabpool;


for i = 1:Ntrials

%gen trajectory
x0 = [x_lb] + [x_ub-x_lb]*rand(Nx,1);
U = u_lb + (u_ub-u_lb)*rand(Nu,len-1);
%with point mass dynamics only:
for k = 1:len
if(k==1)
traj(:,k) = x0;
else
traj(:,k) = A*traj(:,k-1) + B*U(:,k-1);

for dim=1:Nx
if(traj(dim,k)>x_ub)
traj(dim,k) = x_ub;
end %endif
if(traj(dim,k)<x_lb)
traj(dim,k) = x_lb;
end %endif

end %endfor

end %endif
end %endfor
Traj_all(:,:,i) = traj;
U_all(:,:,i) = U;

end %endfor



parfor i = 1:Ntrials
    i

%traj = .25*randn(dim,len);
trajec = Traj_all(:,:,i);
%get robustness going
r_exact = alt_getRobustnessP_vector(trajec,P,wavparams,1);
r_approx = alt_getRobustnessP_vector(trajec,P,wavparams,0);
err(i) = r_exact-r_approx;
err_rel(i) = err(i)/abs(r_exact);
end
[mean((err_rel)) std((err_rel))]
[mean(err) std(err)]
rmse_compute(err)
figure;
hist(err);grid on;
%save('empirical_err/RobTest_N500_AnP.mat','Traj_all','U_all','err_rel','err','len');
