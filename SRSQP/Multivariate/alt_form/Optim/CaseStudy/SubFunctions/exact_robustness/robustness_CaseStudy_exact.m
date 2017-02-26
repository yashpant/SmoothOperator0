function [r_P,r_P_der] = robustness_CaseStudy_exact(X,ExactParams)
% robustness for CaseStudy
%Phi = eventually(p_1 in Tset) AND eventually(p_2 in Tset) AND always( not(p_1 in
% Z_1) OR (z_1 in Rule1)) AND always( not(p_2 in Z_1) OR (z_2 in Rule1)) AND
% always( not(p_1 in Z_2) OR (z_1 in Rule2)) AND always( not(p_2 in Z_2) OR (z_2 in Rule2))
% AND always(not(p_1 in NoFly)) AND always(not(p_2 in NoFly)) AND
% always( norm(p_1-p_2)^2 <= dmin^2)

%trajectory Sx (n*N) (n = sys dim, N = steps)
%%
r_P_der = [];
r_phi = zeros(9,1);

Nx = ExactParams.dim_x;
Nu = ExactParams.dim_u;
N = ExactParams.len;

state_1 = reshape(X(1:Nx*N),Nx,N);
state_2 = reshape(X(Nx*N+(N-1)*Nu+1:Nx*N+(N-1)*Nu+Nx*N),Nx,N);

pos_1 = state_1(4:6,:);
pos_2 = state_2(4:6,:);
%%
r_phi(1) = robustness_eventually_P_exact(pos_1,ExactParams.Terminal);
r_phi(2) = robustness_eventually_P_exact(pos_2,ExactParams.Terminal);
r_phi(7) = robustness_always_notP_exact(pos_1,ExactParams.NoFly);
r_phi(8) = robustness_always_notP_exact(pos_2,ExactParams.NoFly);

r_phi(9) = robustness_always_SafeDist_exact([pos_1;pos_2],ExactParams.d_min);

r_phi(3) = robustness_always_IfAthenB_exact(pos_1,ExactParams.Zone1, ...
    ExactParams.Zone1_rules);
r_phi(4) = robustness_always_IfAthenB_exact(pos_2,ExactParams.Zone1, ...
    ExactParams.Zone1_rules);
r_phi(5) = robustness_always_IfAthenB_exact(pos_1,ExactParams.Zone2, ...
    ExactParams.Zone2_rules);
r_phi(6) = robustness_always_IfAthenB_exact(pos_2,ExactParams.Zone2, ...
    ExactParams.Zone2_rules);


r_P = min(r_phi);