%% load
load('TestData_opt7.mat');
xtemp = x;
load('TestData_opt4.mat');
x1 = x;
x = xtemp;
%% plot
dim = optParams.dim;
P_feas = optParams.P_feas;
P_final = optParams.P_final;
len = optParams.len;
P1 = optParams.P1;
if(dim<=3)
    figure;
    plot(P_feas,'Color','gray','Alpha',0.7);
    hold on;
    plot(P1,'Color','red','Alpha',0.7);
    hold on;
    plot(P_final,'Color','green','Alpha',0.7);
    hold on;
    traj_x = reshape(x(1:dim*len),dim,len);
    traj_x0 = reshape(x_0(1:dim*len),dim,len);
    traj_x1 = reshape(x1(1:dim*len),dim,len);
    hold on;
    
    for i = 1:len
        plot(traj_x0(1,i),traj_x0(2,i),'bo');hold on;
        plot(traj_x(1,i),traj_x(2,i),'k*');hold on;
        plot(traj_x1(1,i),traj_x1(2,i),'rd');hold on;
        if(i==2)
            legend('Feasible set','Unsafe set','Terminal Set','Init. Traj.',...
                'Traj. \gamma=0.1','Traj. \gamma=0.001');
        end
    end
    grid on;
end
xlabel('x_1');
ylabel('x_2');

