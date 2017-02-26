clc;close all;
disp('Initializing problem...');
genCodeCostFn = 0;
nQuads = 2;
%% sys dynamics
g = 9.8;
m = 0.5;
h = 1/5; %5Hz

tempA = zeros(6,6);
%tempA(1:3,4:6) = eye(3);
tempA(4:6,1:3) = eye(3);
tempB = zeros(6,3);
%tempB(4,1) = g;
%tempB(5,2) = -g;
%tempB(6,3) = 1/m;
tempB(1,1) = g;
tempB(2,2) = -g;
tempB(3,3) = 1/m;

sys_c = ss(tempA,tempB,eye(6),zeros(6,3));
sys_d = c2d(sys_c,h);

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

%E_dash = double(getcondvects(dim));

%%
disp('Loading params...')

load('Data_case/ControllerData_fine.mat','wp');
wavparams.NoFly = wp{1};
wavparams.Terminal = wp{2};
wavparams.Zone1 = wp{3};
wavparams.Zone2 = wp{4};

load('Data_case/AltitudeRuleCoefs.mat','rule_wp');
wavparams.Zone1_rules = rule_wp{1};
wavparams.Zone2_rules = rule_wp{2};

%% in order to compute exact robustness
Zone1_rules = Polyhedron('lb',1,'ub',xmax-1);
Zone2_rules = Polyhedron('lb',0,'ub',xmax/2);
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
disp('Initial trajectory...');
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
u1_0 = x1_feas.u(:);
u2_0 = x2_feas.u(:);
%plot init trajs
figure(1);
for i = 1:optParams.len
   hold on;
   plot3(x1_feas.z(4,i),x1_feas.z(5,i),x1_feas.z(6,i),'b*');
   hold on 
   plot3(x2_feas.z(4,i),x2_feas.z(5,i),x2_feas.z(6,i),'ro');
end
grid on;
zlabel('z');

%% dynamics in u

% LB_U = -U_feas.b(1)*ones(1,optParams.len-1);
% UB_U = U_feas.b(2)*ones(1,optParams.len-1);


dim_u = size(optParams.B,2);
dim_x = size(optParams.A,1);
dim = dim_x;
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

A_x0 = [eye(dim);A_x0];
B_U = [zeros(optParams.dim_x,(optParams.len-1)*optParams.dim_u);B_U];
%% for 2 quadrotors
optParams.A_x0_nQuad = kron(eye(nQuads),A_x0);
optParams.B_U_nQuad = kron(eye(nQuads),B_U);

LB = -repmat(optParams.U_feas.b(1:dim_u),nQuads*(optParams.len-1),1);
UB = repmat(optParams.U_feas.b(1:dim_u),nQuads*(optParams.len-1),1);
optParams.x0 = [optParams.x0_1;optParams.x0_2];
u_0 = [u1_0;u2_0];

%% project state constraints onto inputs
disp('Getting sets for the optimization...');
Feas_Set = Polyhedron('lb',[-5*ones(5,1);0.1],'ub',5*ones(6,1)); %feas set
X_lims = Polyhedron('lb',repmat([-5*ones(5,1);0.1],optParams.len*nQuads,1), ...
    'ub',repmat(5*ones(6,1),optParams.len*nQuads,1)); %limits on the nQuads
H1 = X_lims.A*optParams.B_U_nQuad;
g1 = X_lims.b-X_lims.A*optParams.A_x0_nQuad*optParams.x0;

H1_sparse = sparse(size(H1,1),size(H1,2));
[rr,cc] = find(H1);
H1_sparse(rr,cc) = H1(rr,cc);

U_new = Polyhedron('A',H1_sparse,'b',g1);
U_lims = Polyhedron('lb',LB,'ub',UB);
U_intersect = intersect(U_new,U_lims);
%U_intersect.minHRep;
load('CaseData/SetsForU_Case.mat','U_intersect'); %minHRep for [0;-2;2] etc init condition

%% Gen. Code
disp('Generating Code...');
if(genCodeCostFn)
    cfg=coder.config('mex');
    arg_ins = {coder.typeof(u_0),coder.typeof(optParams)};
    codegen -config cfg objfun_case_u -report -args arg_ins -o objfun_case_u_mex
end

%% Optimization
disp('Computing Control...');
tic;
options = optimset('Algorithm','sqp','Display','iter','MaxIter',1000,'TolConSQP',1e-6,...
    'UseParallel','always','MaxFunEval',1000000,'ObjectiveLimit',-eps,'GradObj','off');
%options.TolFun = 10^(-10);
%options.TolCon = 10;
%[x,fval,flag] = ...
[u_opt,fval,exitflag,output] = fmincon(@(u)objfun_case_u_mex(u,optParams),u_0,U_intersect.A,U_intersect.b,[],[],[],[], ...
    [],options);
time_taken_mins = toc/60
%% %%

save('CaseData/Case_20_u_sat_maxz.mat','u_opt','u_0','optParams','ExactParams');


%% plot
Nx = optParams.dim_x;
Nu = optParams.dim_u;
N = optParams.len;

X = optParams.A_x0_nQuad*optParams.x0 + optParams.B_U_nQuad*u_opt; %get states
state_1 = reshape(X(1:Nx*N),Nx,N);
state_2 = reshape(X(Nx*N+1:end),Nx,N);

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


%% plot inputs
input_1 = reshape(u_opt(1:numel(u_opt)/nQuads),Nu,N-1);
input_2 = reshape(u_opt(numel(u_opt)/nQuads+1:end),Nu,N-1); %Nx*N+(N-1)*Nu+Nx*N

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



