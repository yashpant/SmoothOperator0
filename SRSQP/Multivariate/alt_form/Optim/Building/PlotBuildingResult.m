%%
close all;figure;
P = Polyhedron('lb',[10 22],'ub',[19 28]);
hold on;
plot(P,'Color','green','Alpha',0.2);
axis([1 24 15 40]);
load('/home/mlab-retro/Documents/AnytimeCPS2015/MATLAB/Miscellaneous/Multivariate/alt_form/Optim/Building/BldgData/SatTraj.mat')
z = SimBldg(u_opt',optParams);
hold all;
plot(z(4,:),'b','linewidth',3)
hold all
plot(Sys.model_data.X(4,1:24),'--r','linewidth',3)
load('/home/mlab-retro/Documents/AnytimeCPS2015/MATLAB/Miscellaneous/Multivariate/alt_form/Optim/Building/BldgData/MaxTraj.mat')
z = SimBldg(u_opt',optParams);
hold all;
plot(z(4,:),'k','linewidth',3)
hold all
plot(Sys.model_data.X(4,1:24),'--m','linewidth',2)
load('/home/mlab-retro/Documents/AnytimeCPS2015/MATLAB/Miscellaneous/Multivariate/alt_form/Optim/Building/BldgData/Sim_hot_exact.mat')
z = SimBldg(x,optParams);
hold all;
plot(z(4,:),':g','linewidth',2);
load('/home/mlab-retro/Documents/AnytimeCPS2015/MATLAB/Miscellaneous/Multivariate/alt_form/Optim/Building/BldgData/Sim_SA_exact.mat')
z = SimBldg(x,optParams);
hold all;
plot(z(4,:),'-.c','linewidth',2);
legend('Comfort Zone','SR-SQP (B)','BS (B)','SR-SQP (R)','BS (R)','R-SQP','SA');
xlabel('Time Steps (hours)');
ylabel('Zone Temperature (C)');