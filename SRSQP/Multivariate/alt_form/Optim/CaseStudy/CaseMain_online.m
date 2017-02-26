% case study main

%general options
genCode =0;
genCodeCostFn = 1;
getCoeffs = 0;
getRuleCoeffs = 0;

%% sys dynamics
g = 9.8;
m = 0.5;
h = 1/5; %5Hz

tempA = sparse(6,6);
%tempA(1:3,4:6) = eye(3);
tempA(4:6,1:3) = eye(3);
tempB = sparse(6,3);
%tempB(4,1) = g;
%tempB(5,2) = -g;
%tempB(6,3) = 1/m;
tempB(1,1) = g;
tempB(2,2) = -g;
tempB(3,3) = 1/m;

sys_c = ss(tempA,tempB,eye(6),zeros(6,3));
sys_d = c2d(sys_c,h);
%% viz some trajectories for sanity check
x = zeros(6,1);
Tmax = 100;
X = zeros(6,Tmax);
for t = 1:Tmax
    u=[randn(2,1);1];
    x = sys_d.A*x + sys_d.B*u;
    X(:,t) = x;
    
end
plot3(X(4,:),X(5,:),X(6,:),'o');grid on
axis('square');
xlabel('x');ylabel('y');zlabel('z');

%% describe sets and their signed dist approxs
xmin = -7;
xmax = 7;
dx = 0.1;

NoFly = Polyhedron('lb',[-1 -1 0],'ub',[1 1 5]);
Terminal = Polyhedron('lb',[3 3 0],'ub',[4 4 1]);
Feasible = Polyhedron('lb',[xmin*ones(1,2) 0],'ub',xmax*ones(1,3));
LimitSet = Polyhedron('lb',[-5*ones(1,2) 0],'ub',5*ones(1,3));
Zone1 = Polyhedron('lb',[xmin*ones(1,2) 0],'ub',[0 xmax xmax]);
Zone2 = Polyhedron('lb',[0 xmin 0],'ub',xmax*ones(1,3));

wp = cell(4,1);

close all;
figure(1);
plot(intersect(Feasible,LimitSet),'Color','gray','Alpha',0.1);
hold on;
plot(intersect(Zone1,LimitSet),'Color','orange','Alpha',0.2);
hold on;
plot(intersect(Zone2,LimitSet),'Color','black','Alpha',0.2);
hold on;
plot(intersect(NoFly,LimitSet),'Alpha',0.3);
hold on;
plot(intersect(Terminal,LimitSet),'Color','green','Alpha',0.3);
xlabel('x');ylabel('y');

%%
dim = size(Feasible.A,2);
grid_x = xmin:dx:xmax;
dist_array_xy = zeros(numel(grid_x)*ones(1,dim));

optional.k_min = -5;
optional.k_max = +5;
optional.j_min = 0;
optional.j_max = 0; %2 works well in R2
E_dash = double(getcondvects(dim));


if(genCode)
    cfg=coder.config('mex');
        
    %for get coeffs (grid_x,dist_array_xy,dx,j_min,j_max,k_min,k_max,E_dash,0);
    arg_ins = {coder.typeof(grid_x),coder.typeof(dist_array_xy),coder.typeof(dx),...
       coder.typeof(optional.j_min),coder.typeof(optional.j_max),...
       coder.typeof(optional.k_min), coder.typeof(optional.k_max), coder.typeof(E_dash),...
       coder.typeof(0)};
    codegen -config cfg alt_getCoefficientsVector_genable -report -args arg_ins
        
end
%%
optional.filename = [];
if(getCoeffs)
   parfor i=1:4
       if(i==1)
    %optional.savefile = 'NoFly_Params.mat'; 
   wavparams_NoFly = WavSignedDistVector_Rn(NoFly,xmin,xmax,dx,0);%,optional);
  wp{i} = wavparams_NoFly;	
       end
        if(i==2)
   %optional.savefile = 'Terminal_Params.mat';
   wavparams_Terminal = WavSignedDistVector_Rn(Terminal,xmin,xmax,dx,0);%,optional);
   wp{i} = wavparams_Terminal;	
        end
        if(i==3)
        wavparams_Zone1 = WavSignedDistVector_Rn(Zone1,xmin,xmax,dx,0);%,optional);    
	wp{i} = wavparams_Zone1;
        end
        if(i==4)
        wavparams_Zone2 = WavSignedDistVector_Rn(Zone2,xmin,xmax,dx,0);%,optional);    
	wp{i} = wavparams_Zone2;
        end
   end
end
if(getCoeffs)
save('ControllerData_fine.mat','wp');
%save('Controller data','wavparams_NoFly','wavparams_Terminal','wavparams_Zone1',...
%    'wavparams_Zone2');
else
    load('ControllerData_fine.mat','wp');
    wavparams.NoFly = wp{1};
    wavparams.Terminal = wp{2};
    wavparams.Zone1 = wp{3};
    wavparams.Zone2 = wp{4};
end
%% flight (altitude) rules for Zone1 and Zone2
%figure(2);
Zone1_rules = Polyhedron('lb',1,'ub',xmax-1);
Zone2_rules = Polyhedron('lb',0,'ub',xmax/2);
%plot(Zone1_rules,'Color','red');
%hold on;
%plot(Zone2_rules,'Color','blue');
%grid on;

rule_wp = cell(2,1);

if(getRuleCoeffs)
    parfor i = 1:2
    if(i==1)
        rule_wp{i} = WavSignedDistVector_Rn_noGen(Zone1_rules,-2,xmax,dx,1);
    end
    if(i==2)
        rule_wp{i} = WavSignedDistVector_Rn_noGen(Zone2_rules,-2,xmax,dx,1);
    end
    end
    save('Data_case/AltitudeRuleCoefs.mat','rule_wp');
else
    load('Data_case/AltitudeRuleCoefs.mat','rule_wp');
    wavparams.Zone1_rules = rule_wp{1};
    wavparams.Zone2_rules = rule_wp{2};
end
%%
ExactParams.Terminal = Terminal;
ExactParams.NoFly = NoFly;
ExactParams.Zone1 = Zone1;
ExactParams.Zone2 = Zone2;
ExactParams.Zone1_rules = Zone1_rules;
ExactParams.Zone2_rules = Zone2_rules;
ExactParams.dim_x = 6;
ExactParams.dim_u = 3;
ExactParams.len = 20;
ExactParams.d_min = 0.2;
%% optim para
optParams.wavparams = wavparams;
optParams.d_min = 0.2;
optParams.dim_x = 6;
optParams.dim_u = 3;
optParams.len = 20;
optParams.A = sys_d.A;
optParams.B = sys_d.B;

max_ang = deg2rad(30);
max_thr = 1.5;
temp = Polyhedron('lb',[-max_ang -max_ang -max_thr],'ub',[max_ang max_ang max_thr]);
optParams.U_feas.A = temp.A;
optParams.U_feas.b = temp.b;

temp = Polyhedron('lb',[-5*ones(5,1);0.1],'ub',5*ones(6,1)); %feas set
optParams.P_feas.A = temp.A;
optParams.P_feas.b = temp.b;

Terminal_velocities = Polyhedron('lb',-1*ones(3,1),'ub',1*ones(3,1));
temp = Terminal_velocities*Terminal;
optParams.P_final.A = temp.A;
optParams.P_final.b = temp.b; 

optParams.robCost = 1;
optParams.robConstr = 0;
optParams.gamma = 10^(-2);

%% initialize optimization
'Getting feas traj'
x1 = zeros(optParams.len*optParams.dim_x+ (optParams.len-1)*optParams.dim_u,1); %for one quad
x1(4:6) = [-2 2 2];
x1_feas = getFeasTraj_case(x1,optParams);
x2 = zeros(optParams.len*optParams.dim_x+ (optParams.len-1)*optParams.dim_u,1); %for one quad
x2(4:6) = [-2 -2 2];
x2_feas = getFeasTraj_case(x2,optParams);

% init states for the optimization
optParams.x0_1 = x1(1:6);
optParams.x0_2 = x2(1:6);

x_0 =[x1_feas.x0;x2_feas.x0]; %starting point of opt (a trajectory)
%plot init trajs
figure(1);
for i = 1:optParams.len
   hold on;
   plot3(x1_feas.z(4,i),x1_feas.z(5,i),x1_feas.z(6,i),'g*');
   hold on 
   plot3(x2_feas.z(4,i),x2_feas.z(5,i),x2_feas.z(6,i),'ko');
end
grid on;
zlabel('z');
%% generate some code for cost function
if(genCodeCostFn)
     cfg=coder.config('mex');
       
    %for get coeffs (grid_x,dist_array_xy,dx,j_min,j_max,k_min,k_max,E_dash,0);
    arg_ins = {coder.typeof(x_0),coder.typeof(optParams)};
    codegen -config cfg objfun_case -report -args arg_ins -o objfun_case_mex_onlne
    %codegen -config cfg confun_case -report -args arg_ins -o confun_case_mex
end

%% optimization shite
'Starting optimization'
tic;
options = optimset('Algorithm','sqp','Display','iter','MaxIter',1000,'TolConSQP',1e-6,...
    'UseParallel','always','MaxFunEval',1000000,'GradObj','off');
%options.TolFun = 10^(-10);
%options.TolCon = 10;
%[x,fval,flag] = ...
[x,fval,exitflag,output] = fmincon(@(x)objfun_case_mex_onlne(x,optParams),x_0,[],[],[],[],[],[], ...
    @(x)confun_case(x,optParams),options);
toc;
%% plot
Nx = optParams.dim_x;
Nu = optParams.dim_u;
N = optParams.len;
state_1 = reshape(x(1:Nx*N),Nx,N);
state_2 = reshape(x(Nx*N+(N-1)*Nu+1:Nx*N+(N-1)*Nu+Nx*N),Nx,N);

pos_1 = state_1(4:6,:);
pos_2 = state_2(4:6,:);

for i = 1:size(pos_1,2)
    hold on
    plot3(pos_1(1,i),pos_1(2,i),pos_1(3,i),'g*');
    hold on;
    plot3(pos_2(1,i),pos_2(2,i),pos_2(3,i),'ko');
    pause(0.1)
end

figure;
hold on;
for i = 1:size(pos_1,2)
    hold on;
plot(i,norm(pos_1(:,i)-pos_2(:,i)),'r.')
end
hold on;
plot(0:optParams.len,optParams.d_min*ones(1,optParams.len+1),'g--');
grid on;
%%

save('Case_20_len.mat','x','x_0','optParams','ExactParams');
time_taken_mins = toc/60

%%
input_1 = reshape(x(1+Nx*N:Nx*N+(N-1)*Nu),Nu,N-1);
input_2 = reshape(x(Nx*N+(N-1)*Nu+1+Nx*N:end),Nu,N-1); %Nx*N+(N-1)*Nu+Nx*N

figure;
subplot(311)
plot(rad2deg(input_1(1,:)),'b');
hold on;
plot(rad2deg(input_2(1,:)),'r');
grid on;
subplot(312)
plot(rad2deg(input_1(2,:)),'b');
hold on;
plot(rad2deg(input_2(2,:)),'r');
grid on;
subplot(313)
plot((input_1(3,:)),'b');
hold on;
plot((input_2(3,:)),'r');
grid on;



