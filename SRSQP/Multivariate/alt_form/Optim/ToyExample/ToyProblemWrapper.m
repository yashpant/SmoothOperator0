% wrapper
Nruns = 25;

global rand_x0
rand_x0 = 1;
global display_on;
display_on = 0;
global robustness_max

clear SmoothOpt;
global taliro_SmoothRob; %flag
global SmoothOpt; %params
taliro_SmoothRob = 0;

disp(' ')
disp('The specification:')
phi = '[]_[0,20.0]!a /\ <>_[0,20.0]b'
%phi = '[]_[0,2.0]!a' % /\ <>_[0,2.0]b'
%phi = '<>_[0,3.9]b'
preds(1).str = 'a';
disp('Type "help monitor" to see the syntax of MTL formulas')
P_pred_a = Polyhedron('lb',[-1 -1],'ub',[1 1]); %The polyhedron for predicate a
preds(1).A = P_pred_a.A;
preds(1).b = P_pred_a.b;

preds(2).str = 'b';
disp('Type "help monitor" to see the syntax of MTL formulas')
P_pred_b = Polyhedron('lb',[2 2],'ub',[2.5 2.5]);
preds(2).A = P_pred_b.A;
preds(2).b = P_pred_b.b;

%% SR-SQP sat mode
if(1)
robustness_max = 1;

time_per_shot_SRSQP = zeros(Nruns-1,1);
robustness_one_shot_SRSQP = zeros(Nruns-1,1);
for ii = 1:Nruns
    ii
    main_example2_toy_1shot_u_param;
    if(display_on)
    pause
    end
    if(ii>1)
        time_per_shot_SRSQP(ii-1) = time_taken;
        robustness_one_shot_SRSQP(ii-1) = -fval;
    end
end
disp('Satisfaction rate')
sum(robustness_one_shot_SRSQP>=0)/numel(robustness_one_shot_SRSQP)
disp('Mean, max and std. exec time for SQSQP')
[mean(time_per_shot_SRSQP) max(time_per_shot_SRSQP) std(time_per_shot_SRSQP)]
end
sop_200 = time_per_shot_SRSQP.*(robustness_one_shot_SRSQP>=0);
[mean(sop(find(sop_200))) std(sop(find(sop_200)))]


%% SR-SQP max mode
if(1)
robustness_max = 1;
display_on = 0;
time_per_shot_SRSQP = zeros(Nruns-1,1);
robustness_one_shot_SRSQP = zeros(Nruns-1,1);
for ii = 1:Nruns
    ii
    pause(0.2)
    example2_toy_1shot;
    if(display_on)
    pause
    end
    if(ii>1)
        time_per_shot_SRSQP(ii-1) = time_taken;
        robustness_one_shot_SRSQP(ii-1) = -fval;
        robustness_one_shot_SRSQP_actual(ii-1) = dp_taliro(phi,preds,traj_x',[0:.1:(len-1)/10]',[],[],[])
    end
end
disp('Satisfaction rate')
sum(robustness_one_shot_SRSQP>=0)/numel(robustness_one_shot_SRSQP)
disp('Mean, max and std. exec time for SQSQP')
[mean(time_per_shot_SRSQP) max(time_per_shot_SRSQP) std(time_per_shot_SRSQP)]
end
save('Data/MaxModeTimes_50.mat','time_per_shot_SRSQP','robustness_one_shot_SRSQP','robustness_one_shot_SRSQP_actual','ii');
%% BluSTL
if(1)
time_per_shot_BS = zeros(Nruns-1,1);
robustness_one_shot_BS = zeros(Nruns-1,1);
global time_1shot;
%robustness_one_shot_BS = zeros(Nruns-1,1);
for ii = 1:Nruns
    ToyExample2d;
    if(ii>1)
        ii
        time_per_shot_BS(ii-1) = time_1shot;
        SX = Sys.model_data.X;
        ST = Sys.model_data.time;
        robustness_one_shot_BS(ii-1) = dp_taliro(phi,preds,SX',ST',[],[],[]);
        if(display_on) 
            pause;
        end
        
    end
end

disp('Mean, max and std. exec time for BluSTL')
[mean(time_per_shot_BS) max(time_per_shot_BS) std(time_per_shot_BS)]
end
%% save data
if(0)
save('Data/RobMode_20.mat','time_per_shot_SRSQP','time_per_shot_BS','ii');
end
%% SR-SQP Max Mode
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