%% Defining the plant dynamics 
% The toolbox is organized around one main class, called STLClti. An
% STLClti object is primarily a continuous Linear Time Invariant (LTI) 
% system with inputs, outputs and disturbances. Hence, a constructors for this 
% class takes matrices to define such an LTI. We first define and A and B
% matrices for state evolution:

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
Nx = 6;
Nu = 3;
%% BS Begins
disp('Initializing problem')
A = [sys_c.A zeros(Nx,Nx);zeros(Nx,Nx) sys_c.A];
Bu = [sys_c.B zeros(Nx,Nu);zeros(Nx,Nu) sys_c.B];



%%
% Later on, we will use a disturbance signals so we need to define a Bw 
% matrice. This signal will not influence the state dynamics, though, so we
% set Bw to be 0. 
Bw = zeros(2*Nx,1); 

%%
% Next we define the output dynamics, i.e., C, Du and Dw matrices. Here we
% have a single output $y(t) = x1(t)+w(t)$.
C = eye(Nx*2);
Du = zeros(2*Nx,2*Nu);
Dw = zeros(2*Nx,1);

%%
% Now we can call the main constructor of STLC_lti class. 
Sys= STLC_lti(A,Bu,Bw,C,Du,Dw); 

%%
% In the next section, we will define the different settings for the 
% control synthesis experiment. Before that, we define some initial state:
% if(~exist('rand_x0','var'))
% fixed_x0 = 0;
% else
%    fixed_x0 = ~rand_x0; 
% end
% 
% if(fixed_x0);
% x0 = [-2;-2];
% else %in [-2 -1.1]^2
% x0 = -2 + (-1.1+2)*rand(2,1);    
% end
x0 = zeros(12,1);
x0(4:6) = [-2;2;2];
x0(10:12) = [-2;-2;2];
Sys.x0= x0;

%% Defining the controller
% We start by defining the time instants for the whole experiment, the discrete time
% step ts for the controller and the horizon L in number of time steps.  x
Sys.time = 0:2*h:8; 
Sys.ts=h; % sampling time for controller
Sys.L=20;  % horizon is 2\4s in that case

%%
% Next we declare some constraints on control inputs, here, lower and upper
% bounds:
Sys.u_ub = repmat([deg2rad(30);deg2rad(30);1.5],2,1);  % upper bound on u 
Sys.u_lb = -repmat([deg2rad(30);deg2rad(30);1.5],2,1);   % lower bound on u

%%
% Then the following define a signal temporal logic formula to be satisfied
% by the system. Note that times in the temporal operators are continuous,
% not discrete steps. 


%Sys.stl_list = {'(alw_[0,4] not ((-x4(t)<=1) and (x4(t)<=1) and (-x5(t)<=1) and (x5(t)<=1) and (-x6(t)<=-0.001) and (x6(t)<=5))) and (alw_[0,4] not ((-x10(t)<=1) and (x10(t)<=1) and (-x11(t)<=1) and (x11(t)<=1) and (-x12(t)<=-0.001) and (x12(t)<=5))) and (ev_[0,4]((-x4(t)<=-3) and (x4(t)<=4) and (-x5(t)<=-3) and (x5(t)<=4) and (-x6(t)<=-0.001) and (x6(t)<=1))) and (ev_[0,4]((-x10(t)<=-3) and (x10(t)<=4) and (-x11(t)<=-3) and (x11(t)<=4) and (-x12(t)<=-0.001) and (x12(t)<=1)))'};

%Sys.stl_list = {'(alw_[0,4] not ((-x4(t)<=1) and (x4(t)<=1) and (-x5(t)<=1) and (x5(t)<=1) and (-x6(t)<=-0.001) and (x6(t)<=5))) and (alw_[0,4] not ((-x10(t)<=1) and (x10(t)<=1) and (-x11(t)<=1) and (x11(t)<=1) and (-x12(t)<=-0.001) and (x12(t)<=5))) and (ev_[0,4]((-x4(t)<=-3) and (x4(t)<=4) and (-x5(t)<=-3) and (x5(t)<=4) and (-x6(t)<=-0.001) and (x6(t)<=1))) and (ev_[0,4]((-x10(t)<=-3) and (x10(t)<=4) and (-x11(t)<=-3) and (x11(t)<=4) and (-x12(t)<=-0.001) and (x12(t)<=1))) and (alw_[0,4] not ((abs(x4(t)-x10(t))<=0.2) and (abs(x5(t)-x11(t))<=0.2) and (abs(x6(t)-x12(t))<=0.2))) and (alw_[0,4] (not ((-x4(t)<=7) and (x4(t)<=0) and (-x5(t)<=7) and (x5(t)<=7) and (-x6(t)<=0) and (x6(t)<=5)) or ((-x6(t)<=-1) and (x6(t)<=4))))'};

%Sys.stl_list = {'(alw_[0,4] not ((-x4(t)<=1) and (x4(t)<=1) and (-x5(t)<=1) and (x5(t)<=1) and (-x6(t)<=-0.001) and (x6(t)<=5))) and (alw_[0,4] not ((-x10(t)<=1) and (x10(t)<=1) and (-x11(t)<=1) and (x11(t)<=1) and (-x12(t)<=-0.001) and (x12(t)<=5))) and (ev_[0,4]((-x4(t)<=-3) and (x4(t)<=4) and (-x5(t)<=-3) and (x5(t)<=4) and (-x6(t)<=-0.001) and (x6(t)<=1))) and (ev_[0,4]((-x10(t)<=-3) and (x10(t)<=4) and (-x11(t)<=-3) and (x11(t)<=4) and (-x12(t)<=-0.001) and (x12(t)<=1))) and (alw_[0,4] not ((abs(x4(t)-x10(t))<=0.2) and (abs(x5(t)-x11(t))<=0.2) and (abs(x6(t)-x12(t))<=0.2))) and (alw_[0,4] (not ((-x4(t)<=7) and (x4(t)<=0) and (-x5(t)<=7) and (x5(t)<=7) and (-x6(t)<=0) and (x6(t)<=5)) or ((-x6(t)<=-1) and (x6(t)<=4)))) and (alw_[0,4] (not ((-x10(t)<=7) and (x10(t)<=0) and (-x11(t)<=7) and (x11(t)<=7) and (-x12(t)<=0) and (x12(t)<=5)) or ((-x12(t)<=-1) and (x12(t)<=4))))'};  

%full 2-quad spec.
Sys.stl_list = {'(alw_[0,4] not ((-x4(t)<=1) and (x4(t)<=1) and (-x5(t)<=1) and (x5(t)<=1) and (-x6(t)<=-0.001) and (x6(t)<=5))) and (alw_[0,4] not ((-x10(t)<=1) and (x10(t)<=1) and (-x11(t)<=1) and (x11(t)<=1) and (-x12(t)<=-0.001) and (x12(t)<=5))) and (ev_[0,4]((-x4(t)<=-3) and (x4(t)<=4) and (-x5(t)<=-3) and (x5(t)<=4) and (-x6(t)<=-0.001) and (x6(t)<=1))) and (ev_[0,4]((-x10(t)<=-3) and (x10(t)<=4) and (-x11(t)<=-3) and (x11(t)<=4) and (-x12(t)<=-0.001) and (x12(t)<=1))) and (alw_[0,4] not ((abs(x4(t)-x10(t))<=0.2) and (abs(x5(t)-x11(t))<=0.2) and (abs(x6(t)-x12(t))<=0.2))) and (alw_[0,4] (not ((-x4(t)<=7) and (x4(t)<=0) and (-x5(t)<=7) and (x5(t)<=7) and (-x6(t)<=0) and (x6(t)<=5)) or ((-x6(t)<=-1) and (x6(t)<=4)))) and (alw_[0,4] (not ((-x10(t)<=7) and (x10(t)<=0) and (-x11(t)<=7) and (x11(t)<=7) and (-x12(t)<=0) and (x12(t)<=5)) or ((-x12(t)<=-1) and (x12(t)<=4)))) and (alw_[0,4] (not ((-x4(t)<=0) and (x4(t)<=7) and (-x5(t)<=7) and (x5(t)<=7) and (-x6(t)<=0) and (x6(t)<=5)) or ((-x6(t)<=0) and (x6(t)<=3.5)))) and (alw_[0,4] (not ((-x10(t)<=0) and (x10(t)<=7) and (-x11(t)<=7) and (x11(t)<=7) and (-x12(t)<=0) and (x12(t)<=5)) or ((-x12(t)<=0) and (x12(t)<=3.5)))) '};  
%% 1 shot
if(1)
% Now we are ready to compile the controller for our problem. 

disp('Getting controller')
controller = get_controller(Sys,'boolean')

%%
% Note that by default, the objective function will minimize the 1-norm of
% the input. 

%% Testing the controller
% The simplest mode to run our system with the newly created controller is in
% open loop. This is done with the following command:
disp('Open loop mode')
run_open_loop(Sys, controller);

%%
if(exist('display_on','var'))
    if(display_on)
        disp_plot=1;
    else
        disp_plot=0;close all;
    end
else
    disp_plot = 1;
end
if(disp_plot)
    %%
close all;
xmin = -7;xmax=7;
NoFly = Polyhedron('lb',[-1 -1 0],'ub',[1 1 5]);
Terminal = Polyhedron('lb',[3 3 0],'ub',[4 4 1]);
Feasible = Polyhedron('lb',[xmin*ones(1,2) 0],'ub',xmax*ones(1,3));
LimitSet = Polyhedron('lb',[-5*ones(1,2) 0],'ub',5*ones(1,3));
Zone1 = Polyhedron('lb',[xmin*ones(1,2) 0],'ub',[0 xmax xmax]);
Zone2 = Polyhedron('lb',[0 xmin 0],'ub',xmax*ones(1,3));

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
xlabel('x');ylabel('y');zlabel('z');
hold on;
plot3(Sys.model_data.X(4,:),Sys.model_data.X(5,:),Sys.model_data.X(6,:),'o');
hold on;
plot3(Sys.model_data.X(10,:),Sys.model_data.X(11,:),Sys.model_data.X(12,:),'o');

end
end
%%

if(0) %the MPC part, comment it back in if needed
    
% We can run the system in closed loop, but this is not very interesting,
% because w is 0 anyway. Let us change this: 
Sys.Wref = Sys.time*0.;
%Sys.Wref(30:40) = 1; 
%Sys.Wref(60:80) = -0.5; 


%%
% and the specifications:
%Sys.stl_list = {'alw (ev_[0,1.] alw_[0,0.5] ( abs(y1(t)-w1(t)) < 0.1))'};
%Sys.stl_list = {'alw_[0,4] (abs(y1(t)) < 1)'};
disp('getting MPC');
controller = get_controller(Sys,'boolean');

%%
% This time we will only plot input and outputs, i.e., disable the state
% plotting:
disp('MPC mode')
close;
Sys.plot_x =[]; % default was Sys.plot_x = [1 2]
run_deterministic(Sys, controller);

%%
% More examples are given in the folder BluSTL/examples. In particular the 
% hvac_room case study demonstrates the adversarial scenario, as well as
% plot customization. The idea is to create a class derived from STLC_lti
% and specialize the update_plot method. 
figure;
%plot(P_feas,'Color','gray','Alpha',0.7); % not a constraint in this code
hold on;
plot(P1,'Color','red','Alpha',0.7); %unsafe set
hold on; 
plot(P_final,'Color','green','Alpha',0.7); %terminal set
hold on;
plot(Sys.x0(1),Sys.x0(2),'ko');
hold on;
plot(Sys.model_data.X(1,:),Sys.model_data.X(2,:),'*');grid on
axis([-10 10 -10 10]);
legend('Unsafe set','Terminal set','Initial state','trajectory')

end