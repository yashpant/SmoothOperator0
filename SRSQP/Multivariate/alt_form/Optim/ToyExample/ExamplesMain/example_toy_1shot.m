% example 4, optimization
% (always not P1) AND/OR (min/max) (always not P2)

%% get params
P1 = Polyhedron('lb',[-1 -1],'ub',[1 1]);
xmin = -5;
xmax = 5;
dx = 0.1;

optional_P1.filename = 'P1_j0_P.mat';
optional_P1.savefile = 'P1_j0_P.mat';

[Params_P1] = WavSignedDistVector(P1,xmin,xmax,dx,0,optional_P1);

%pause;
%close all;clc;

%% optimization data
dim = 2;
len = 30;

P_feas = Polyhedron('lb',[-2.5 -2.5],'ub',[2.5 2.5]);
P_final = Polyhedron('lb',[2.0 2.0],'ub',[2.5 2.5]);
P_term = P_final;
%U_feas = Polyhedron('lb',[-0.32 -0.32],'ub',[0.32 0.32]);
%U_feas = Polyhedron('lb',[-2.52 -2.52],'ub',[2.52 2.52]);
U_feas = Polyhedron('lb',[-0.32 -0.32],'ub',[0.32 0.32]);


optParams.P_final = P_final;
optParams.U_feas = U_feas;
%system
optParams.A = eye(2);
optParams.B = eye(2);
% robustness in obj and/or constr
optParams.robCost = 1;
optParams.robConstr = 0;

x0 = [-2;-2];
%x0 = [-1.5;0];
optParams.gamma = 10^(-2);
optParams.x0 = x0;
optParams.dim = dim;
optParams.len = len;
optParams.Params_P1 = Params_P1;
optParams.P1 = P1;
optParams.P_feas = P_feas;
%% start opt
clc;
% init traj gen
if(1) %via reg MPC
    x_0 = [x0;rand((len-1)*dim,1);rand((len-1)*size(optParams.B,2),1)];
    x_feas = getFeasTraj(x_0,optParams);
    x_0 = x_feas.x0;
    if(sum(isnan(x_0))>0)
        'x_0 infeasible'
        keyboard;
    end
else
    load('TestData_opt4','x') %from somewhere else
    x_0 = x;
end
%x_0 = [x0;[-1.5;1];[-0.5
'got init traj'
%%  optim
tic;
options = optimset('Algorithm','sqp','Display','iter','MaxIter',1000,'TolConSQP',1e-6,...
    'UseParallel',true,'MaxFunEval',1000000,'GradObj','off');
%options.TolFun = 10^(-10);
%options.TolCon = 10;
%[x,fval,flag] = ...
[x,fval,exitflag,output] = fmincon(@(x)objfun_toy(x,optParams),x_0,[],[],[],[],[],[], ...
    @(x)confun_toy(x,optParams),options);


save('Data/TestData_opt10.mat','x','x_0','optParams');
time_taken_mins = toc/60

%% plot
dim = optParams.dim;
P_feas = optParams.P_feas;
P_final = optParams.P_final;
len = optParams.len;
P1 = optParams.P1;
if(dim<=3)
    figure;
    plot(P_feas,'Color','gray','Alpha',0.7);
    hold on;
    plot(P1,'Color','red','Alpha',0.7);
    hold on;
    plot(P_final,'Color','green','Alpha',0.7);
    hold on;
    traj_x = reshape(x(1:dim*len),dim,len);
    traj_x0 = reshape(x_0(1:dim*len),dim,len);
    hold on;
    
    for i = 1:len
        plot(traj_x0(1,i),traj_x0(2,i),'bo');hold on;
        plot(traj_x(1,i),traj_x(2,i),'k*');hold on;
        
        if(i==2)
            legend('Feasible set','Unsafe set','Terminal Set','Init. Traj.',...
                'Traj. \gamma=0.1');
        end
    end
    grid on;
end


