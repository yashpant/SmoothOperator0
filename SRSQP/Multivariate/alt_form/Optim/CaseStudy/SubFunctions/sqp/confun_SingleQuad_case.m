function [c, ceq] = confun_SingleQuad_case(x,x0,optParams)
%%
%x = [z;u]
%z+ = Az + bu

Nx = optParams.dim_x;
N = optParams.len;
Nu = size(optParams.B,2);

% load('OptParams.mat');

% robustness
rho = 0;
if(optParams.robConstr)
%traj = reshape(x(1:Nx*N),Nx,N);
%rho = alt_getRobustnessP_vector(traj,optParams.P1,optParams.Params_P1,0);
rho = 0;
end


% initial condition
ceq_init  = x(1:6,1) - x0(1:6,1);

% dynamics;
bigA_zz = [zeros(Nx,Nx*N); kron(eye(N-1),optParams.A) zeros(Nx*N-Nx,Nx)];
bigA_zu = [zeros(Nx,(N-1)*Nu); kron(eye(N-1),optParams.B)];
bigA_uu = eye((N-1)*Nu);zeros((N-1)*Nu);
bigA_uz = zeros((N-1)*Nu, N*Nx);
bigA = [bigA_zz bigA_zu;bigA_uz bigA_uu];
ceq_dym = x-bigA*x;
ceq_dym = ceq_dym(Nx+1:end);

% state feasibility
A_big_x = kron(eye(N),optParams.P_feas.A);
B_big_x = repmat(optParams.P_feas.b,N,1);
c_state = A_big_x*x(1:Nx*N) - B_big_x; 

% input feasibility
A_big_u = kron(eye(N-1),optParams.U_feas.A);
B_big_u = repmat(optParams.U_feas.b,N-1,1);
c_input = A_big_u*x(Nx*N+1:end) - B_big_u;

% terminal set
%c_term = optParams.P_final.A*x(N*Nx-Nx+1:N*Nx)-optParams.P_final.b;
c_term = 0; %if no term set here
if(optParams.robConstr)
c = [c_state;c_input;c_term;-rho];
else
c = [c_state;c_input;c_term]; 
end

ceq = [ceq_init;ceq_dym];