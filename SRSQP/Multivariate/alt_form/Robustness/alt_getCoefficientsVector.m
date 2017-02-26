function [C_00k,D_ejk] = alt_getCoefficientsVector(grid_x,dist_array_xy,dx,j_min,j_max,k_min,k_max,E_dash,viz)

%see Nonlinear Approximation, R. A. DeVore, Acta Numerica 1998, Cambridge
%University Press
%%
E = E_dash(2:end,:);
dim = log2(size(E,1)+1);

K = zeros(dim^(numel(k_min:k_max)),dim);

for i = 1:dim^(numel(k_min:k_max))
   K(i,:) = permn2(k_min:k_max,dim,i); 
end
%K = permn(k_min:k_max,dim);
k_sz = size(K,1);

X = zeros(dim^numel(grid_x),dim);

for i = 1:dim^(numel(grid_x))
   X(i,:) = permn2(grid_x,dim,i); 
end
%X = permn(grid_x,dim);

x_sz = size(X,1);
j_sz = numel(j_min:j_max);
D_ejk = zeros(size(E_dash,1),j_sz,k_sz);
C_00k = zeros(1,1,k_sz);

phis = zeros(1,dim);
psis = zeros(1,dim);
%%
%compute C_00k
for k_dims = 1:k_sz
    summation_C = 0;
    for x_dims = 1:x_sz
        x = (2^0)*X(x_dims,:)' - K(k_dims,:)';
        
        %[phis,~] = arrayfun(@(t) MeyerWavelet(t),x);
        %phis = (2^0)*phis';
        
        
        for dimen = 1:dim
           phis(dimen) = (2^0)*MeyerWavelet(x(dim));
        end
        
        bigPhi_00k_x = prod(phis);
        g_x = dist_array_xy(x_dims);
        summation_C = summation_C + (dx^dim)*g_x*bigPhi_00k_x; %check use of dim and dx here
    end
    C_00k(1,1,k_dims) = summation_C;
end

for e = 1:size(E,1)
    j_ix = 0;
    for j=j_min:1:j_max
        j_ix = j_ix + 1;
        for k_dims = 1:k_sz %k along each dim
            summation_D = 0;
            for x_dims = 1:x_sz
                if(viz)
                [e j k_dims x_dims];
                end
                x = (2^j)*X(x_dims,:)' - K(k_dims,:)';
                
                for dimen = 1:dim
                   [phis(dimen),psis(dimen)] = MeyerWavelet(x(dimen));
                   
                end
                
                %[phis,psis] = arrayfun(@(t) MeyerWavelet(t),x);
                phis = (2^j)*phis;
                psis = (2^j)*psis;
                %keyboard
                phi_e = (1-E(e,:)).*phis;
                psi_e = E(e,:).*psis;
                bigPsi_temp = phi_e+psi_e;
                bigPsi_ejk_x = prod(bigPsi_temp);
                % check if vectorization preserves locations in this %
                g_x = dist_array_xy(x_dims);
                summation_D = summation_D + (dx^dim)*g_x*bigPsi_ejk_x; %check use of dim and dx here
            end
            D_ejk(e,j_ix,k_dims) = summation_D;
        end
    end
    
end



