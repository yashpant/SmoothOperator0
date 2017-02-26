%%
%load('Wavparams2_j0_k10_ToyExample.mat');
dim = optParams.dim;
P_feas = AuxParams.P_feas;
P_final = AuxParams.P_final;
len = optParams.len;
P_unsafe = AuxParams.P_unsafe;
if(dim<=3)
    figure;
    plot(P_feas,'Color','gray','Alpha',0.7);
    hold on;
    plot(P_unsafe,'Color','red','Alpha',0.7);
    hold on;
    plot(P_final,'Color','green','Alpha',0.7);
    hold on;
    grid on;
end

%% for eventually;
x_pt = [1;-1]
plot(x_pt(1),x_pt(2),'*');
[f2,g2] = robustness_eventually_P(x_pt,optParams.Params_P_term,1);
sd_true = SignedDist(x_pt,optParams.P_final.A,optParams.P_final.b);
[f2 sd_true]
[g2 [robustness_eventually_P_genable_parallel(x_pt+[1;0],optParams.Params_P_term,1)-sd_true;robustness_eventually_P_genable_parallel(x_pt+[0;1],optParams.Params_P_term,1)-sd_true] ... 
    [SignedDist(x_pt+[1;0],optParams.P_final.A,optParams.P_final.b)-sd_true;SignedDist(x_pt+[0;1],optParams.P_final.A,optParams.P_final.b)-sd_true]]


%% for eventually;
clc;
x_pt = [1;1]
offset = 0.1;
plot(x_pt(1),x_pt(2),'*');
[f1,g1] = alt_getRobustnessP_and_der_vector_genable(x_pt,optParams.Params_P_unsafe,1);
f1_true = SignedDist(x_pt,optParams.P_unsafe.A,optParams.P_unsafe.b);
[f1 f1_true]
[g1 [alt_getRobustnessP_and_der_vector_genable(x_pt+offset*[1;0],optParams.Params_P_unsafe,1)-f1_true;alt_getRobustnessP_and_der_vector_genable(x_pt+offset*[0;1],optParams.Params_P_unsafe,1)-f1_true] ... 
    [SignedDist(x_pt+offset*[1;0],optParams.P_unsafe.A,optParams.P_unsafe.b)-f1_true;SignedDist(x_pt+offset*[0;1],optParams.P_unsafe.A,optParams.P_unsafe.b)-f1_true]]

%%
[f,g] = objfun2_toy_using_mex(x,optParams);
figure;
stairs(g,'*');
offset = 1;
g_finite = zeros(numel(g),1);
for i = 1:numel(g)
    e = zeros(size(g_finite));
    e(i) = offset;
    g_finite(i) = (objfun2_toy_using_mex(x+e,optParams)-f)/norm(e);
    
end
hold all;
stairs(g_finite,'.-');
%%
[f,g] = main_objfun2_u_toy_using_mex(u_0,optParams);
figure;
stairs(g);
offset = 0.1;
g_finite = zeros(numel(u_0),1);
for i = 1:optParams.dim_u*(optParams.len-1)
    e = zeros(size(g_finite));
    e(i) = offset;
    g_finite(i) = (main_objfun2_u_toy_using_mex(u_0+e,optParams) - f)/norm(e);
    
end
hold all;
stairs(g_finite);

%%
figure;
[f,g] = objfun2_toy_using_mex(x,optParams);
stairs(g);
hold all
u = x(41:end);
[fu,gu] = main_objfun2_u_toy_using_mex(u,optParams);
stairs(gu);
figure;
stairs(optParams.B_U'*g(1:40)-gu);

%% check states

state(:,1) = optParams.x0;
uu = reshape(u_opt,2,19);
for i = 2:20
    state(:,i) = optParams.A*state(:,i-1) + optParams.B*uu(:,i-1);
end
figure;
plot(state(1,:),state(2,:),'o');hold all;
state_fast = optParams.A_x0*optParams.x0 + optParams.B_U*u_opt;
state_fast_rs = reshape(state_fast,2,20);
grid on;plot(state_fast_rs(1,:),state_fast_rs(2,:),'k*');
