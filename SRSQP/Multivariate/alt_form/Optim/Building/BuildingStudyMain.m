%Bldg main
%%
Ts = 60*60; %sampling time, 1 hr
H = 24; %24 hour, 1 step per hour
[A,B,Bd,C] = linear_model(Ts); 
disturbances = sim_disturbance(1999,4,1,1999,4,2); %disturbances for a day
U_feas = Polyhedron('lb',-500,'ub',2000);
P_feas = Polyhedron('lb',0*ones(4,1),'ub',50*ones(4,1));
dec_factor = 10;
P1_comfort = Polyhedron('lb',22,'ub',28);
I1=10:19; %hours of interest in the 24 hour period
optParams.I1 = I1;
load('ComfortParams.mat');
optParams.wavparams = wavparams;
%% generate some code
genCodeCostFn = 0;
if(genCodeCostFn)
    cfg=coder.config('mex');
    x = 0;
    %for get coeffs (grid_x,dist_array_xy,dx,j_min,j_max,k_min,k_max,E_dash,0);
    arg_ins = {coder.typeof(x),coder.typeof(wavparams.C),coder.typeof(wavparams.D),...
         coder.typeof(wavparams.k_min),coder.typeof(wavparams.k_max), ...
         coder.typeof(wavparams.j_min),coder.typeof(wavparams.j_max)};
    codegen -config cfg getWavApprox -report -args arg_ins -o getWavApprox_mex
    %codegen -config cfg confun_case -report -args arg_ins -o confun_case_mex
end

%% init optimization
x_0 = 21*ones(4,1); %init state
optParams.x0 = x_0;
optParams.disturbances = disturbances;
optParams.dim = 4;
optParams.len = 24;
optParams.dim_d = 3;
optParams.dim_u = 1;
optParams.A = A;
optParams.B = B*10; %check this
optParams.Bd = Bd;
optParams.P_feas.A = P_feas.A;
optParams.P_feas.b = P_feas.b;
optParams.U_feas.A = U_feas.A;
optParams.U_feas.b = U_feas.b;
optParams.P1_comfort.A = P1_comfort.A;
optParams.P1_comfort.b = P1_comfort.b;
optParams.exact = 0;
%% get feas traj
x_feas = bldg_getFeasTraj(x_0,optParams,disturbances);
u_0 = x_feas.u;
%% gen code for objfun and confun

if(1)
    disp('Code gen');
    cfg=coder.config('mex');
    arg_ins = {coder.typeof(u_0),coder.typeof(optParams),coder.typeof(I1)};
    codegen -config cfg objfun_bldg_gen -report -args arg_ins
    %objfun_bldg_gen(u,optParams,I1)
end


%% optimization

if(1)

'Starting optimization'
tic;
options = optimset('Algorithm','sqp','Display','iter','MaxIter',1000,'TolConSQP',1e-6,...
    'UseParallel','always','MaxFunEval',1000000,'GradObj','off');
%options.TolFun = 10^(-10);
%options.TolCon = 10;
%[x,fval,flag] = ...

%[x,fval,exitflag,output] = fmincon(@(x)objfun_bldg_gen_mex(x,optParams,I1),u_0,[],[],[],[],[],[], ...
%    @(x)confun_bldg(x,optParams),options);

[x,fval,exitflag,output] = fmincon(@(x)objfun_bldg_gen(x,optParams,I1),u_0,[],[],[],[], ...
    -500*ones(size(u_0)),2000*ones(size(u_0)),[],options);


toc

else %simann
optParams.exact = 1;
tic
[x,fmin,Nruns,info]=sa_mincon_bldg(0.8,optParams,u_0);
end
%%
plot_bldg(x,optParams);
%%
save('BldgData/Sim_SA_exact.mat','x','optParams');

