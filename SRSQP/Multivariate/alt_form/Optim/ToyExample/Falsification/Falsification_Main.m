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
len = 20;

P_feas = Polyhedron('lb',[-2.5 -2.5],'ub',[2.5 2.5]);
P_final = Polyhedron('lb',[2.0 2.0],'ub',[2.5 2.5]);
P_term = P_final;
P_init = Polyhedron('lb',[-2.5 -2.5],'ub',[-1.25 2.5]);
optParams.P_init = P_init;
%U_feas = Polyhedron('lb',[-0.32 -0.32],'ub',[0.32 0.32]);
%U_feas = Polyhedron('lb',[-2.52 -2.52],'ub',[2.52 2.52]);
U_feas = Polyhedron('lb',[-0.52 -0.52],'ub',[0.52 0.52]);

optParams.P_final = P_final;
optParams.U_feas = U_feas;
%system
A = eye(2);B = eye(2);
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

%% plot sets
figure(1);
plot(P_feas,'Color','gray','Alpha',0.3);
hold on;
plot(P_term,'Color','green','Alpha',0.7);
hold on;
plot(P1,'Color','red','Alpha',0.7);
hold on;
plot(P_init,'Color','blue','Alpha',0.7);

%% dynamics and some testing
[K,S,E] = dlqr(A,B,eye(2),eye(2),zeros(2,2)); %optimal, aggressive
[K,prec] = place(A,B,[.7;.7]); %less aggressive
x0 = [-2.25;+2.25];
x = zeros(2,1:len);
x(:,1) = x0;
x_star = [2.25;2.25];

optParams.x0 = x0;
optParams.x_star = x_star;
optParams.K = K;
% change of fucking vars
% z = x-x_star
% v = u-u_star

u_star = inv(B)\(eye(2)-A)*x_star;
%u_star = 0;
% z+ = Az+Bv = (A-BK)z
% i.e. x+-x_star = (A-BK)*(x-x_star)
% i.e. x+ = (A-BK)(x-x_star) + x_star 
% since A=I, x+ = (A-B*K)*x + B*K*x_star

for t = 2:len
    x(:,t) = (A-B*K)*x(:,t-1) + B*K*x_star;%(A-B*K)*(x(:,t-1) - x_star) + x_star;
       
end

hold on;
plot(x(1,:),x(2,:),'bo');

% making sure of this
optParams.Params_P1.j_max = 0;

%% gen code if necessary
genCodeCostFn = 0;
if(genCodeCostFn)
    cfg=coder.config('mex');
    Params_P1.j_max = 0;  
    %for get coeffs (grid_x,dist_array_xy,dx,j_min,j_max,k_min,k_max,E_dash,0);
    arg_ins = {coder.typeof(x),coder.typeof(Params_P1)};
    codegen -config cfg robustness_toy_genable -report -args arg_ins -o robustness_toy_mex
    %codegen -config cfg confun_case -report -args arg_ins -o confun_case_mex
end

%% optimization
tic;
options = optimset('Algorithm','sqp','Display','iter','MaxIter',1000,'TolConSQP',1e-6,...
    'UseParallel','always','MaxFunEval',1000000,'GradObj','off');
exact = 1; %use what function
[x0_opt,fval,exitflag,output] = fmincon(@(x)objfun_toy_falsification(x,optParams,exact),x0,[],[],[],[],[],[], ...
    @(x)confun_toy_falsification(x,optParams),options);


%save('Data/TestData_opt10.mat','x','x_0','optParams');
%time_taken_mins = toc/60

%% plot

%sets


% traj
x = zeros(2,len);
x(:,1) = x0_opt;
x(:,1) = x_sr;
for t = 2:len
    x(:,t) = (A-B*K)*x(:,t-1) + B*K*x_star;%(A-B*K)*(x(:,t-1) - x_star) + x_star;
       
end

hold on;
plot(x(1,:),x(2,:),'kd');

%% legend shite 
h  =legend('Feasible set','Terminal set','Unsafe set','X_0','\mathbb{x}_0','\mathbb{x}_{\tilde{\rho}}','7');
set(h,'Interpreter','latex');
h  =legend('Feasible set','Terminal set','Unsafe set','X_0','Init. pt. (and traj.)','Init. state (and traj.) SQP ', ...
    'Init. state (and traj.) SA') ;
xlabel('x_1');ylabel('x_2');

