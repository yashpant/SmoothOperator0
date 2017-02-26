%load('Case_20_u_nofevlim.mat'); %z(5,:)
%load('Case_20_u_zmax.mat'); %-z(6,:)
%load('Case_20_u_0.mat'); %-z(4,:)
%% plot sets
xmin = -7;
xmax = 7;
dx = 0.1;

NoFly = Polyhedron('lb',[-1 -1 0],'ub',[1 1 5]);
Terminal = Polyhedron('lb',[3 3 0],'ub',[4 4 1]);
Feasible = Polyhedron('lb',[xmin*ones(1,2) 0],'ub',xmax*ones(1,3));
LimitSet = Polyhedron('lb',[-5*ones(1,2) 0],'ub',5*ones(1,3));
Zone1 = Polyhedron('lb',[xmin*ones(1,2) 0],'ub',[0 xmax xmax]);
Zone2 = Polyhedron('lb',[0 xmin 0],'ub',xmax*ones(1,3));



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

%% plot trajs
Nx = optParams.dim_x;
Nu = optParams.dim_u;
N = optParams.len;

%1
load('Case_20_u_nofevlim.mat'); %z(5,:)
%load('Case_20_u_zmax.mat'); %-z(6,:)
%load('Case_20_u_0.mat'); %-z(4,:)


X = optParams.A_x0_nQuad*optParams.x0 + optParams.B_U_nQuad*u_opt; %get states
state_1 = reshape(X(1:Nx*N),Nx,N);
state_2 = reshape(X(Nx*N+1:end),Nx,N);

pos_1 = state_1(4:6,:);
pos_2 = state_2(4:6,:);

figure(1)
hold all;
plot3(pos_1(1,:),pos_1(2,:),pos_1(3,:),'--*b','linewidth',2,'MarkerSize',8)
hold all;
plot3(pos_2(1,:),pos_2(2,:),pos_2(3,:),'--og','linewidth',2,'MarkerSize',8)

%2

%load('Case_20_u_nofevlim.mat'); %z(5,:)
load('Case_20_u_zmax.mat'); %-z(6,:)
%load('Case_20_u_0.mat'); %-z(4,:)


X = optParams.A_x0_nQuad*optParams.x0 + optParams.B_U_nQuad*u_opt; %get states
state_1 = reshape(X(1:Nx*N),Nx,N);
state_2 = reshape(X(Nx*N+1:end),Nx,N);

pos_1 = state_1(4:6,:);
pos_2 = state_2(4:6,:);


figure(1)
hold all;
plot3(pos_1(1,:),pos_1(2,:),pos_1(3,:),'--*r','linewidth',2,'MarkerSize',8)
hold all;
plot3(pos_2(1,:),pos_2(2,:),pos_2(3,:),'--oc','linewidth',2,'MarkerSize',8)

%3

%load('Case_20_u_nofevlim.mat'); %z(5,:)
%load('Case_20_u_zmax.mat'); %-z(6,:)
load('Case_20_u_0.mat'); %-z(4,:)


X = optParams.A_x0_nQuad*optParams.x0 + optParams.B_U_nQuad*u_opt; %get states
state_1 = reshape(X(1:Nx*N),Nx,N);
state_2 = reshape(X(Nx*N+1:end),Nx,N);

pos_1 = state_1(4:6,:);
pos_2 = state_2(4:6,:);


figure(1)
hold all;
plot3(pos_1(1,:),pos_1(2,:),pos_1(3,:),'--*m','linewidth',2,'MarkerSize',8)
hold all;
plot3(pos_2(1,:),pos_2(2,:),pos_2(3,:),'--ok','linewidth',2,'MarkerSize',8)

%
set(gca,'Xtick',-5:5:5)
set(gca,'Ytick',-5:5:5)
set(gca,'Ztick',0:5:5)
legend('Airspace','Zone1','Zone2','NoFly','Terminal','p_{1}^1','p_{2}^1','p_{1}^2','p_{2}^2','p_{1}^3','p_{2}^3')

%%
figure;
hold on;
for i = 1:size(pos_1,2)
    hold on;
plot(i,norm(pos_1(:,i)-pos_2(:,i)),'r.')
end
hold on;
plot(0:optParams.len,optParams.d_min*ones(1,optParams.len+1),'g--');
grid on;

%% get rob
X = [state_1(:);u_opt(1:numel(u_opt)/2);state_2(:);u_opt(numel(u_opt)/2+1:end)];
robustness_CaseStudy_exact(X,ExactParams)
-objfun_case_u_mex(u_opt,optParams)-optParams.d_min
