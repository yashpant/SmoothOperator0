% 2d example
% (always not Unsfafe) and (eventually Terminal)
%% get params
codegeneration = 0;
no_disturbances = 0;

clc;close all;
disp('Initializing problem...');


% prelim set stuff
U_feas = Polyhedron('lb',-500,'ub',2000);
xmin = 0;
xmax = 50;
P_feas = Polyhedron('lb',xmin*ones(4,1),'ub',xmax*ones(4,1));
dec_factor = 10;
P1_comfort = Polyhedron('lb',22,'ub',28);
optParams.P1_comfort.A = P1_comfort.A;
optParams.P1_comfort.b = P1_comfort.b;
clear SmoothOpt;
global taliro_SmoothRob; %flag
global SmoothOpt; %params
taliro_SmoothRob = 0;
genCode = 0;
preds(1).A = P1_comfort.A;
preds(1).b = P1_comfort.b;

load('ComfortParams.mat');
SmoothOpt.preds.Wavparams = wavparams;
optParams.wavparams = wavparams;


%% optimization data
disp('Control problem data...');

Ts = 60*60; %sampling time, 1 hr
len = 24; %24 hour, 1 step per hour
[A,B,Bd,C,~] = linear_model(Ts); 
disturbances = sim_disturbance(1999,4,1,1999,4,2); %disturbances for a day

if(no_disturbances)
   disturbances.d1 = zeros(size(disturbances.d1));
   disturbances.d2 = zeros(size(disturbances.d2));
   disturbances.d3 = zeros(size(disturbances.d3));
end
optParams.disturbances = disturbances;
Nd = size(Bd,2);
D = [];
for i = 1:len-1
    D((i-1)*Nd+1:i*Nd,1) = [disturbances.d1(i);disturbances.d2(i);disturbances.d3(i)];
end
optParams.D = D; %vector of disturbances;




dim = 4;
dim_u = 1;
dim_d = 3;
optParams.dim = dim;
optParams.len = len;
optParams.dim_u = dim_u;
optParams.dim_d = dim_d;

I1=10:19; %hours of interest in the 24 hour period
optParams.I1 = I1;


AuxParams.U_feas = U_feas;
optParams.U_feas.A = U_feas.A;
optParams.U_feas.b = U_feas.b;
% system dynamics
optParams.A = A;
optParams.Bd = Bd;
optParams.B = B*10; %check this
% robustness in obj and/or constr
optParams.robCost = 1;
optParams.robConstr = 0;

%% init state

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

optParams.P_feas.A = P_feas.A;
optParams.P_feas.b = P_feas.b;

%% redo constraints for control input

LB_U = -U_feas.b(1)*ones(1,optParams.len-1);
UB_U = U_feas.b(2)*ones(1,optParams.len-1);
dim_u = size(optParams.B,2);
% translate state constraints if needed
A_x0 = [];
for i = 1:optParams.len-1
    A_x0 = [A_x0;optParams.A^i];    
end
%
clear B_U
for i = 1:optParams.len-1
    for j = 1:optParams.len-1
        if(i>=j)
            B_U((i-1)*dim+1:i*dim,(j-1)*dim_u+1:j*dim_u) = optParams.A^(i-j)*optParams.B;
        else
            B_U((i-1)*dim+1:i*dim,(j-1)*dim_u+1:j*dim_u) = zeros(size(optParams.A*optParams.B));     
        end
       %B_U((i-1)*dim+1:i*dim,(j-1)*dim_u+1:j*dim_u) = ((i-j)>=0)*optParams.A^(i-j)*optParams.B + ((i-j)<0)*zeros(size(optParams.A*optParams.B));     
    end
end

clear B_D
for i = 1:optParams.len-1
    for j = 1:optParams.len-1
        if(i>=j)
            B_D((i-1)*dim+1:i*dim,(j-1)*dim_d+1:j*dim_d) = optParams.A^(i-j)*optParams.Bd;
        else
            B_D((i-1)*dim+1:i*dim,(j-1)*dim_d+1:j*dim_d) = zeros(size(optParams.A*optParams.Bd));
        end
       %B_D((i-1)*dim+1:i*dim,(j-1)*dim_d+1:j*dim_d) = ((i-j)>=0)*optParams.A^(i-j)*optParams.Bd + ((i-j)<0)*zeros(size(optParams.A*optParams.Bd));     
    end
end

% optParams.A_x0 = A_x0;
% optParams.B_U = B_U;
optParams.A_x0 = [eye(dim);A_x0];
optParams.B_U = [zeros(optParams.dim,(optParams.len-1)*optParams.dim_u);B_U];
optParams.B_D = [zeros(optParams.dim,(optParams.len-1)*optParams.dim_d);B_D]; % all distbs stacked up

%% start opt
disp('Getting init trajectory...');
% init traj gen
if(1) %via reg MPC
   x_feas = bldg_getFeasTraj(x0,optParams,disturbances);
   u_0 = x_feas.u';
    if(sum(isnan(u_0))>0)
        'x_0 infeasible'
        keyboard;
    end
else
    load('TestData_opt4','x') %from somewhere else
    x_0 = x;
end
close all;

%%
disp('Mapping all constraints to inputs...');
X_lims = Polyhedron('lb',(xmin/2)*ones(optParams.len*optParams.dim,1),'ub',(xmax/2)*ones(optParams.len*optParams.dim,1));
H1 = X_lims.A*optParams.B_U;
g1 = X_lims.b-X_lims.A*optParams.A_x0*optParams.x0-X_lims.A*optParams.B_D*D;
U_new = Polyhedron('A',H1,'b',g1);
U_lims = Polyhedron('lb',LB_U,'ub',UB_U);
U_intersect = intersect(U_new,U_lims);
U_intersect.minHRep;



%% gen code for objfun and confun

if(codegeneration)
    disp('Code gen...');
    cfg=coder.config('mex');
    arg_ins = {coder.typeof(u_0),coder.typeof(optParams),coder.typeof(I1),coder.typeof(1),coder.typeof(zeros(optParams.dim,optParams.len))};
    codegen -config cfg objfun_bldg_new_rh -report -args arg_ins
    %objfun_bldg_gen(u,optParams,I1)
end

%%  optim
if(exist('robustness_max','var'))
    if(~isempty(robustness_max))
        objLim = -10*robustness_max - eps*(1-robustness_max);
    else
        objLim = -eps;
    end
    
    else
    objLim = -10;
end

%%
clear options time_taken u_applied dists;
disp('Robustness maximization...')

temp_dists = optParams.B_D*D;
temp_dists = 0.05*reshape(temp_dists,optParams.dim,optParams.len); %5 percent of dists

options = optimset('Algorithm','sqp','Display','off','MaxIter',1000,'TolConSQP',1e-2,...
    'ObjectiveLimit',objLim,'UseParallel','always','MaxFunEval',1000000,'GradObj','on'); %rep 'always' by true
Aeq = zeros(numel(u_0));
Beq = zeros(numel(u_0),1);
x_upto_t = zeros(optParams.dim,optParams.len);
x_t = optParams.x0;
for t=1:optParams.len-1
    t
    x_upto_t(:,t) = x_t; 
    if(t>1)
       u_0 = u_opt;
       Beq(t-1) = u_applied(t-1);
       Aeq(t-1,t-1) = 1;
    end    
    tic;
[u_opt,fval,exitflag,output] = fmincon(@(u)objfun_bldg_new_rh_mex(u,optParams,I1,t,x_upto_t),u_0,[],[],Aeq,Beq,LB_U',UB_U',[],options);
time_taken(t) = toc;
u_applied(t) = u_opt(t);
dists(:,t) = (-temp_dists(:,t)/2) + temp_dists(:,t).*rand(optParams.dim,1);
x_t = optParams.A*x_t + optParams.B*u_applied(t) + dists(:,t);





end
x_upto_t(:,end) = x_t; %last step
%% plot static
clear xx
for i = 1:optParams.len
    if(i==1)
   xx(:,i) = x0; 
    else
   xx(:,i) = optParams.A*xx(:,i-1) + optParams.B*u_applied(i-1) + dists(:,i-1);     
    end     
end

plot(xx(4,:));grid on;
%%
disp('robustness=')
min([28-xx(4,10:19) xx(4,10:19)-22])
