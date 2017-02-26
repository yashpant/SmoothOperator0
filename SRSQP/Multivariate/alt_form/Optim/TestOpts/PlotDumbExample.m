%% plot dumb shite
load('sqp_dumbexample.mat');
figure;

hold on;
plot(x_0(1),x_0(2),'r+','MarkerSize',16,'LineWidth',2);
hold on;
plot(x(1),x(2),'kp','MarkerSize',16,'LineWidth',2);
hold on;
plot(2.5,2.5,'g*','MarkerSize',16,'LineWidth',2);

for i = 1:size(x_samples,2)
   hold on;
   plot(x_samples(1,i),x_samples(2,i),'o','MarkerSize',10);
end
P_feas = Polyhedron('lb',[-2.5 -2.5],'ub',[2.5 2.5]);
P_set = Polyhedron('lb',[-1 -1],'ub',[1 1]);
plot(P_feas,'Color','gray','Alpha',0.3);hold on;
plot(P_set,'Color','red','Alpha',0.3);hold on;
%%
legend('x_0','x','x^*','samples')
xlabel('x_1');
ylabel('x_2');
