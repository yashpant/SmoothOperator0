function plot_bldg(u,optParams)

distbs = [optParams.disturbances.d1'; optParams.disturbances.d2'; ...
    optParams.disturbances.d3'];

z = SimBldg(u,optParams);

figure(4)
hold all;
plot(z(4,:),'linewidth',2);
grid on;
%%
 P = Polyhedron('lb',[10 22],'ub',[19 28]);
 hold on;
 plot(P,'Color','green','Alpha',0.6);
 %%
 legend('Initial','Method','SQP','SA','Comfort temp. (and occupancy time)')
%%
figure(1)
for i = 411:414
    subplot(i)
plot(z(i-410,:));
grid on;
end
grid on;
figure(2)
for i = 411:414
    subplot(i)
    if(i==411)
        plot(u);
        grid on;
    else
        j = i-411;
        stairs(distbs(j,:));
        grid on;
    end
end
grid on;

%%
   I1 = optParams.I1; 
   figure(1);
   subplot(414)
   hold on;
   plot(optParams.I1,-1*optParams.P1_comfort.b(1)*ones(numel(I1),1),'r-');
   hold on
   plot(optParams.I1,1*optParams.P1_comfort.b(2)*ones(numel(I1),1),'r-');
%% 
if(0) %for ws shit
    figure;
    z = SimBldg(u_opt_sat',optParams);
    plot([0:23],z(4,:),'linewidth',2);
    hold all;
    z = SimBldg(u_opt_max',optParams);
    plot([0:23],z(4,:),'linewidth',2);
    hold all;
    plot([0:23],Sys_sat_BS.system_data.X(4,1:24),'linewidth',2);
    hold all;
    plot([0:23],Sys_rmax_BS.system_data.X(4,1:24),'linewidth',2);
     P = Polyhedron('lb',[9 22],'ub',[18 28]);
     hold on;
 plot(P,'Color','green','Alpha',0.6);

xlabel('Time (24 hour format)');
ylabel('Room Temp. (C)')
legend('SRSQP sat','SRSQP rob. max','BS sat','BS rob. max','Comfort temp during 9am-6pm')

end