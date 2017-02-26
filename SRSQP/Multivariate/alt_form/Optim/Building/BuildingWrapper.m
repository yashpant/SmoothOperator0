% wrapper
Nruns = 51;

global rand_x0
rand_x0 = 1;
global display_on;
display_on = 0;
global robustness_max

clear SmoothOpt;
global taliro_SmoothRob; %flag
global SmoothOpt; %params
taliro_SmoothRob = 0;

Ts = 3600; %s
Horizon = 24 %hours
T_index = [0:Ts:(Horizon-1)*Ts];
T_occupied = [9*Ts:Ts:18*Ts] %working hours
disp(' ')
disp('The specification:')
%phi = '[]_[0,2.0]!a /\ <>_[0,2.0]b'
%phi = '[]_[0,2.0]!a' % /\ <>_[0,2.0]b'
phi = '[]_[32400,64800]a'
preds(1).str = 'a';
disp('Type "help monitor" to see the syntax of MTL formulas')
P_pred_a = Polyhedron('lb',22,'ub',28); %The polyhedron for predicate a
preds(1).A = P_pred_a.A;
preds(1).b = P_pred_a.b;

%preds(2).str = 'b';
%disp('Type "help monitor" to see the syntax of MTL formulas')
%P_pred_b = Polyhedron('lb',[2 2],'ub',[2.5 2.5]);
%preds(2).A = P_pred_a.A;
%preds(2).b = P_pred_a.b;

%% SR-SQP sat mode
if(1)
robustness_max = 1;

time_per_shot_SRSQP = zeros(Nruns-1,1);
robustness_one_shot_SRSQP = zeros(Nruns-1,1);
robustness_one_shot_SRSQP_actual = zeros(Nruns-1,1);
for ii = 1:Nruns
    BuildingMain_u_variable;
    if(display_on)
    pause
    end
    if(ii>1)
        time_per_shot_SRSQP(ii-1) = time_taken;
        robustness_one_shot_SRSQP(ii-1) = -fval;
        z = SimBldg(u_opt',optParams);
        robustness_one_shot_SRSQP_actual(ii-1) = dp_taliro(phi,preds,z(4,:)',T_index',[],[],[]);
    end
end
%%
disp('Satisfaction rate')
sum(robustness_one_shot_SRSQP>=0)/numel(robustness_one_shot_SRSQP)
disp('Mean, max and std. exec time for SQSQP')
[mean(time_per_shot_SRSQP) max(time_per_shot_SRSQP) std(time_per_shot_SRSQP)]
disp('95 percentile exec. time')
prctile(time_per_shot_SRSQP,.95)
disp('resulting robustness stats')
[mean(robustness_one_shot_SRSQP) std(robustness_one_shot_SRSQP)]
end
%% BluSTL
if(1)
time_per_shot_BS = zeros(Nruns-1,1);
robustness_one_shot_BS = zeros(Nruns-1,1);
global time_1shot;
%robustness_one_shot_BS = zeros(Nruns-1,1);
for ii = 1:Nruns
    BldgExample;
    if(ii>1)
        ii
        time_per_shot_BS(ii-1) = time_1shot;
        SX = Sys.model_data.X(4,:);
        ST = Sys.model_data.time;
        robustness_one_shot_BS(ii-1) = dp_taliro(phi,preds,SX',ST',[],[],[]);
        if(display_on) 
            pause;
        end
        
    end
end
%%
disp('Mean, max and std. exec time for BluSTL')
[mean(time_per_shot_BS) max(time_per_shot_BS) std(time_per_shot_BS)]
disp('95 percentile exec. time')
prctile(time_per_shot_BS,.95)
disp('resulting robustness stats')
[mean(robustness_one_shot_BS) std(robustness_one_shot_BS)]
end
%% save data
if(1)
save('BldgData/MaxModeTimes.mat','time_per_shot_SRSQP','time_per_shot_BS','ii','robustness_one_shot_BS','robustness_one_shot_SRSQP','robustness_one_shot_SRSQP_actual');
end

%% SA
for ii = 1:Nruns

if(~exist('rand_x0','var'))
fixed_x0 = 1;
else
   fixed_x0 = ~rand_x0; 
end

if(fixed_x0);
x0 = 21*ones(4,1); %init state
else %in [-2 -1.1]^2
x0 = 20 + (22-20)*rand(dim,1);    
end
%x0 = [-1.5;0];
optParams.gamma = 0;%10^(-2);
optParams.x0 = x0;
optParams.exact = 1;
optParams.P_feas.A = P_feas.A;
optParams.P_feas.b = P_feas.b;

disp('Getting init trajectory');
% init traj gen

   x_feas = bldg_getFeasTraj(x0,optParams,disturbances);
   u_0 = x_feas.u;
    if(sum(isnan(u_0))>0)
        'x_0 infeasible'
        keyboard;
    end
close all;


    %[bestsol,fmin,N,info]=sa_mincon_bldg(alpha,optParams,u_0)
tic;[u_opt,fmin,~,~] = sa_mincon_bldg(0.8,optParams,u_0);
time_sa_bldg(ii) = toc;
 z = SimBldg(u_opt,optParams);
         robustness_one_shot_SA_actual(ii) = -fmin;
end
%%
if(0)
%% plotting
load('Data/TimingSat.mat');
L = 20:10:50;
figure;
errorbar(L,[t_20_sr t_30_sr t_40_sr t_50_sr], ...
    [er_20_sr er_30_sr er_40_sr er_50_sr]);
hold all;grid on;
errorbar(L,[t_20_bs t_30_bs t_40_bs t_50_bs], ...
    [er_20_bs er_30_bs er_40_bs er_50_bs]);
xlabel('Horizon Length');
ylabel('Exec. time one shot (s)');
legend('SR-SQP','BluSTL');
title('Sat. always safe and eventually terminal');


%% plotting
L = 20:10:50;
fname_prefix = 'Data/SatModeTimes_eventually_';
for jj = 1:numel(L)
    fname_full = strcat(fname_prefix,num2str(L(jj)));
    load(fname_full);
    mean_BS(jj) = mean(time_per_shot_BS);
    er_BS(jj) = std(time_per_shot_BS);
    mean_SRSQP(jj) = mean(time_per_shot_SRSQP);
    er_SRSQP(jj) = std(time_per_shot_SRSQP);
end
figure;
errorbar(L,mean_SRSQP,er_SRSQP);hold all;
grid on;
errorbar(L,mean_BS,er_BS);hold all;
xlabel('Horizon Length');
ylabel('Exec. time one shot (s)');
legend('SR-SQP','BluSTL');
title('Sat. eventually terminal set');


%% plotting always not
L = 20:10:50;
fname_prefix = 'Data/SatModeTimes_alwaysNot_';
for jj = 1:numel(L)
    fname_full = strcat(fname_prefix,num2str(L(jj)));
    load(fname_full);
    mean_BS(jj) = mean(time_per_shot_BS);
    er_BS(jj) = std(time_per_shot_BS);
    mean_SRSQP(jj) = mean(time_per_shot_SRSQP);
    er_SRSQP(jj) = std(time_per_shot_SRSQP);
end
figure;
errorbar(L,mean_SRSQP,er_SRSQP);hold all;
grid on;
errorbar(L,mean_BS,er_BS);hold all;
xlabel('Horizon Length');
ylabel('Exec. time one shot (s)');
legend('SR-SQP','BluSTL');
title('Sat. always not unsafe set');
end