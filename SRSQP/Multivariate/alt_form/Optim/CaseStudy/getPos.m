function [pos_1,pos_2] = getPos(x,optParams)
Nx = optParams.dim_x;
Nu = optParams.dim_u;
N = optParams.len;
state_1 = reshape(x(1:Nx*N),Nx,N);
state_2 = reshape(x(Nx*N+(N-1)*Nu+1:Nx*N+(N-1)*Nu+Nx*N),Nx,N);

pos_1 = state_1(4:6,:);
pos_2 = state_2(4:6,:);