% 2d example
% (always not Unsfafe) and (eventually Terminal)
%% get params

%clc;close all;
close all;
%clear all;close all;

disp('Initializing problem');
P_unsafe = Polyhedron('lb',[-1 -1],'ub',[1 1]);
P_term = Polyhedron('lb',[2 2],'ub',[2.5 2.5]);

xmin = -5;
xmax = 5;
dx = 0.25;

% prelim set stuff
clear SmoothOpt;
global taliro_SmoothRob; %flag
global SmoothOpt; %params
taliro_SmoothRob = 0;
genCode = 0;
preds(1).A = P_unsafe.A;
preds(1).b = P_unsafe.b;
preds(2).A = P_term.A;
preds(2).b = P_term.b;

for i = 1:numel(preds)
    SmoothOpt.preds.Sets(i).A = preds(i).A;
    SmoothOpt.preds.Sets(i).b = preds(i).b;
    SmoothOpt.preds.WavParams(i).k_min = -20;
    SmoothOpt.preds.WavParams(i).k_max = 20;
    SmoothOpt.preds.WavParams(i).j_min = 0;
    SmoothOpt.preds.WavParams(i).j_max = 0;
    SmoothOpt.preds.WavParams(i).x_min = -5;
    SmoothOpt.preds.WavParams(i).x_max = 5;
    SmoothOpt.preds.WavParams(i).dx = 0.25;
    %SmoothOpt.preds.WavParams(i) = wavparams;
end

%% get params
disp('Wavelet params');
getParams = 0;
if(getParams) %do this once, disable flag and load params next time
wavparams_all = getWaveletParameters(SmoothOpt,genCode);    
    for i = 1:numel(preds)
        
        SmoothOpt.preds.WavParams(i).E = wavparams_all(i).E;
        SmoothOpt.preds.WavParams(i).E_dash = wavparams_all(i).E_dash;
        SmoothOpt.preds.WavParams(i).D_ejk = wavparams_all(i).D_ejk;
        SmoothOpt.preds.WavParams(i).C_00k = wavparams_all(i).C_00k;
        
    end
    save('Wavparams2_j0_k20_ToyExample.mat','SmoothOpt','wavparams_all');
else
    load('Wavparams2_j0_k5_ToyExample.mat');
end


%pause;
%close all;clc;

%% optimization data
disp('Control problem data');
dim = 2;
len = 200;
dim_u = 2;
optParams.dim = dim;
optParams.len = len;
optParams.dim_u = dim_u;

u_lim = 0.52;

P_feas = Polyhedron('lb',[-2.5 -2.5],'ub',[2.5 2.5]);
U_feas = Polyhedron('lb',-u_lim*ones(1,dim),'ub',u_lim*ones(1,dim));

AuxParams.P_final = P_term;
optParams.P_final.A = P_term.A;
optParams.P_final.b = P_term.b;

AuxParams.U_feas = U_feas;
optParams.U_feas.A = U_feas.A;
optParams.U_feas.b = U_feas.b;
% system dynamics
optParams.A = eye(2);
optParams.B = eye(2);
% robustness in obj and/or constr
optParams.robCost = 1;
optParams.robConstr = 0;

%% init state
if(~exist('rand_x0','var'))
fixed_x0 = 0;
else
   fixed_x0 = ~rand_x0; 
end

if(fixed_x0);
x0 = [-2;-2];
else %in [-2 -1.1]^2
x0 = -2 + (-1.1+2)*rand(2,1);    
end
%x0 = [-1.5;0];
optParams.gamma = 0;%10^(-2);
optParams.x0 = x0;
optParams.Params_P_unsafe = SmoothOpt.preds.WavParams(1);
optParams.Params_P_term = SmoothOpt.preds.WavParams(2);

optParams.P_feas.A = P_feas.A;
optParams.P_feas.b = P_feas.b;
optParams.P_unsafe.A = P_unsafe.A;
optParams.P_unsafe.b = P_unsafe.b;
AuxParams.P_unsafe = P_unsafe;
AuxParams.P_feas = P_feas;

%% redo constraints for control input
% H_U = kron(eye(optParams.len-1),optParams.U_feas.A);
% g_U = repmat(optParams.U_feas.b,optParams.len-1,1);
%sanity check
% stairs(H_U*u_0<=g_U);
%or
LB_U = repmat(-u_lim*ones(size(optParams.B,2),1),optParams.len-1,1);
UB_U = repmat(u_lim*ones(size(optParams.B,2),1),optParams.len-1,1);
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
       B_U((i-1)*dim+1:i*dim,(j-1)*dim_u+1:j*dim_u) = (i-j>=0)*optParams.A^(i-j)*optParams.B + (i-j<0)*zeros(size(optParams.A*optParams.B));     
    end
end
% optParams.A_x0 = A_x0;
% optParams.B_U = B_U;
optParams.A_x0 = [eye(dim);A_x0];
optParams.B_U = [zeros(optParams.dim,(optParams.len-1)*optParams.dim_u);B_U];
    
X_lims = Polyhedron('lb',(xmin/2)*ones(optParams.len*optParams.dim,1),'ub',(xmax/2)*ones(optParams.len*optParams.dim,1));
H1 = X_lims.A*optParams.B_U;
g1 = X_lims.b-X_lims.A*optParams.A_x0*optParams.x0;
U_new = Polyhedron('A',H1,'b',g1);
U_lims = Polyhedron('lb',LB_U,'ub',UB_U);
U_intersect = intersect(U_new,U_lims);
%U_intersect.minHRep;

%% start opt
disp('Getting init trajectory');
% init traj gen
if(1) %via reg MPC
    x_0 = [x0;rand((len-1)*dim,1);rand((len-1)*size(optParams.B,2),1)];
    x_feas = getFeasTraj(x_0,optParams);
    x_0 = x_feas.x0;
    u_0 = x_feas.u(:);
    if(sum(isnan(x_0))>0)
        'x_0 infeasible'
        keyboard;
    end
else
    load('TestData_opt4','x') %from somewhere else
    x_0 = x;
end

%% random initialization
rd_u0 = 0;
if(1)
   u_0 = -0.25+(.5)*rand(optParams.dim_u*(optParams.len-1),1); %random u_0; 
   rd_u0 = 1;
   x_0 = [optParams.A_x0*optParams.x0 + optParams.B_U*u_0;u_0];
end




%% gen code for objfun and confun

if(0)
    disp('Code gen');
    CodeGeneratorForOptim;
end
%% sanity cheque
if(0)
tic;objfun2_toy_using_mex(x_0,optParams);toc
tic;objfun2_u_toy_using_mex(u_0,optParams);toc
end
%keyboard
%%  optim
if(exist('robustness_max','var'))
    if(~isempty(robustness_max))
        objLim = -10*robustness_max - eps*(1-robustness_max);
    else
        objLim = -eps;
    end
    
    else
    objLim = -0.24;
end

%%
clear options;
disp('Robustness maximization')
tic;
options = optimset('Algorithm','sqp','Display','iter','MaxIter',1000,'TolConSQP',1e-2,...
    'ObjectiveLimit',objLim,'UseParallel','always','MaxFunEval',1000000,'GradObj','on'); %rep 'always' by true
%options.TolFun = 10^(-10);
%options.TolCon = 10;
%[x,fval,flag] = ...
% [u_opt,fval,exitflag,output] = fmincon(@(u)main_objfun2_u_toy_using_mex(u,optParams),u_0,[],[],[],[],LB_U,UB_U,[],options);
[u_opt,fval,exitflag,output] = fmincon(@(u)main_objfun2_u_toy_using_mex(u,optParams),u_0,U_intersect.A,U_intersect.b,[],[],[],[],[],options);

save('TestData_toyexample2_u_grad.mat','u_opt','u_0','optParams','AuxParams','SmoothOpt');
time_taken = toc;

%% plot
if(exist('display_on','var'))
    if(display_on)
        disp_plot=1;
    else
        disp_plot=0;
    end
else
    disp_plot = 1;
end

if(disp_plot)
    
x_init_and_onwards = optParams.A_x0*optParams.x0 + optParams.B_U*u_opt;
x = [x_init_and_onwards;u_opt];
if(rd_u0) %random u
    x_0 = [optParams.A_x0*optParams.x0 + optParams.B_U*u_0;u_0];
end

dim = optParams.dim;

P_feas = AuxParams.P_feas;
P_final = AuxParams.P_final;
len = optParams.len;
P_unsafe = AuxParams.P_unsafe;
if(dim<=3)
    figure;
%     plot(P_feas,'Color','gray','Alpha',0.7);
%     hold on;
    plot(P_unsafe,'Color','red','Alpha',0.7);
    hold on;
    plot(P_final,'Color','green','Alpha',0.7);
    hold on;
    traj_x = reshape(x(1:dim*len),dim,len);
    traj_x0 = reshape(x_0(1:dim*len),dim,len);
    hold on;
    
    for i = 1:len
        plot(traj_x0(1,i),traj_x0(2,i),'bo');hold on;
        plot(traj_x(1,i),traj_x(2,i),'k*');hold on;
        pause(.1);
        if(i==2)
            legend('Unsafe set','Terminal Set','Init. Traj.',...
                'Traj. \gamma=0.1');
        end
    end
    grid on;
end
axis([-5 5 -5 5]);

end