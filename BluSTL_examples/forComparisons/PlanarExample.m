%% Defining the plant dynamics 
% The toolbox is organized around one main class, called STLClti. An
% STLClti object is primarily a continuous Linear Time Invariant (LTI) 
% system with inputs, outputs and disturbances. Hence, a constructors for this 
% class takes matrices to define such an LTI. We first define and A and B
% matrices for state evolution:

disp('Initializing problem')
A = [1 0 ;
     0 1];
Bu = [1 0;0 1]; %this is the DT system we want

temp.sysd=ss(A,Bu,eye(2),zeros(2,2),0.1);
temp.sysc = d2c(temp.sysd); %get the corresponding CT system

A = temp.sysc.a; %now use the CT matrices for the rest 
Bu = temp.sysc.b;
clear temp;


%%
% Later on, we will use a disturbance signals so we need to define a Bw 
% matrice. This signal will not influence the state dynamics, though, so we
% set Bw to be 0. 
Bw = zeros(2,1); 

%%
% Next we define the output dynamics, i.e., C, Du and Dw matrices. Here we
% have a single output $y(t) = x1(t)+w(t)$.
C = eye(2);
Du = zeros(2,2);
Dw = zeros(2,1);

%%
% Now we can call the main constructor of STLC_lti class. 
Sys= STLC_lti(A,Bu,Bw,C,Du,Dw); 

%%
% In the next section, we will define the different settings for the 
% control synthesis experiment. Before that, we define some initial state:

x0 = [-2;-2];
Sys.x0= x0;

%% Defining the controller
% We start by defining the time instants for the whole experiment, the discrete time
% step ts for the controller and the horizon L in number of time steps.  x
Sys.time = 0:.2:10; 
Sys.ts=.1; % sampling time for controller
Sys.L=50;  % horizon is 2s in that case

%%
% Next we declare some constraints on control inputs, here, lower and upper
% bounds:
Sys.u_ub = [0.52 0.52];  % upper bound on u 
Sys.u_lb = [-0.52 -0.52]; % lower bound on u

%%
% Then the following define a signal temporal logic formula to be satisfied
% by the system. Note that times in the temporal operators are continuous,
% not discrete steps. 

Sys.stl_list = {'(alw_[0,5] not ((abs(y1(t)) <= 1) and (abs(y2(t)) <= 1))) and (ev_[0,5]((-y1(t)<=-2) and (y1(t)<=2.5) and (-y2(t)<=-2) and (y2(t)<=2.5)))'}; %always safe and eventually terminal set

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

%% plot trajectory and sets

close all;
P_feas = Polyhedron('lb',[-2.5 -2.5],'ub',[2.5 2.5]);
P_final = Polyhedron('lb',[2.0 2.0],'ub',[2.5 2.5]);
P_term = P_final;
P1 = Polyhedron('lb',[-1 -1],'ub',[1 1]);
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