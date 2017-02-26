%% Defining the plant dynamics 
% The toolbox is organized around one main class, called STLClti. An
% STLClti object is primarily a continuous Linear Time Invariant (LTI) 
% system with inputs, outputs and disturbances. Hence, a constructors for this 
% class takes matrices to define such an LTI. We first define and A and B
% matrices for state evolution:


% please comment out Sys=Sys.update_plot(); and drawnow; in the STLC_run_open_loop.m file

disp('Initializing problem...')
addpath(genpath('model')); %add the hamlab folder to path
Ts = 60*60; %sampling time, 1 hr
len = 24; %24 hour, 1 step per hour
[A,Bu,Bd,C,sysc] = linear_model(Ts);  %this is the DT system we want
Bu = Bu*10; %get some more control authority
disturbances = sim_disturbance(1999,4,1,1999,4,3); %disturbances, not used



temp.sysd=ss(A,Bu,eye(size(A,1)),zeros(size(A,1),size(Bu,2)),Ts);
temp.sysc = d2c(temp.sysd); %get the corresponding CT system

A = temp.sysc.a; %now use the CT matrices for the rest 
Bu = temp.sysc.b;
clear temp;


%%
% Later on, we will use a disturbance signals so we need to define a Bw 
% matrice. This signal will not influence the state dynamics, though, so we
% set Bw to be 0. 

%temp.sysd=ss(A,Bd,eye(size(A,1)),zeros(size(A,1),size(Bd,2)),Ts);
%temp.sysc = d2c(temp.sysd); %get the corresponding CT system

Bw = zeros(size(A,1),1); 
%Bw = sysc.b(:,[1,3:4]);
%Sys.Wref = [disturbances.d1 disturbances.d2 disturbances.d3];

%%
% Next we define the output dynamics, i.e., C, Du and Dw matrices. Here we
% have a single output $y(t) = x1(t)+w(t)$.
C = eye(size(A,1));
Du = zeros(size(A,1),1);
Dw = zeros(size(A,1),1);

%%
% Now we can call the main constructor of STLC_lti class. 
Sys= STLC_lti(A,Bu,Bw,C,Du,Dw); 

%%
% In the next section, we will define the different settings for the 
% control synthesis experiment. Before that, we define some initial state:

%initial state, random or fixed
if(~exist('rand_x0','var'))
fixed_x0 = 0;
else
   fixed_x0 = ~rand_x0; 
end

if(fixed_x0);
x0 = 21*ones(4,1); %init state
else %in [-2 -1.1]^2
x0 = 20 + (22-20)*rand(size(A,1),1);    
end

Sys.x0= x0;

%% Defining the controller
% We start by defining the time instants for the whole experiment, the discrete time
% step ts for the controller and the horizon L in number of time steps.  x
Sys.time = 0:Ts:(len*2*Ts); 
Sys.ts=Ts; % sampling time for controller
Sys.L=len;  % horizon is 2s in that case

%%
% Next we declare some constraints on control inputs, here, lower and upper
% bounds:
Sys.u_ub = 2000;  % upper bound on u 
Sys.u_lb = -500; % lower bound on u

%%
% Then the following define a signal temporal logic formula to be satisfied
% by the system. Note that times in the temporal operators are continuous,
% not discrete steps. 
Sys.stl_list = {'(alw_[32400,64800]((-x4(t)<-22) and (x4(t)<28)))'};% temp always between 22 and 28C in these times

%% 1 shot
if(1)
% Now we are ready to compile the controller for our problem. 

disp('Getting controller')
controller = get_controller(Sys,'robust')

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
T_start = find(Sys.model_data.time==32400);
T_end = find(Sys.model_data.time==64800);

plot(Sys.model_data.time/Ts,Sys.model_data.X(4,:));%,Sys.model_data.X(2,:),'*');grid on
hold on;
plot(Sys.model_data.time(T_start:T_end)/Ts,22*ones(size(T_start:T_end)),'.-');
hold on;
plot(Sys.model_data.time(T_start:T_end)/Ts,28*ones(size(T_start:T_end)),'.-');
grid on
axis([0 48 0 50]);
xlabel('Time 9hours)');
ylabel('room temp');
legend('Room temp','upper lim','lower lim');
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
