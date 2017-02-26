function C_ejk = getCoefficientsVector(grid_x,dist_array_xy,dx,j_min,j_max,k_min,k_max,E,viz)

%see Nonlinear Approximation, R. A. DeVore, Acta Numerica 1998, Cambridge
%University Press
%%
dim = log2(size(E,1)+0);
K = permn(k_min:k_max,dim);
k_sz = size(K,1);
X = permn(grid_x,dim);
x_sz = size(X,1);
j_sz = numel(j_min:j_max);
C_ejk = zeros(size(E,1),j_sz,k_sz);

%%
for e = 1:size(E,1)
    j_ix = 0;
    for j=j_min:1:j_max
        j_ix = j_ix + 1;
        for k_dims = 1:k_sz %k along each dim
            summation_C = 0;
            parfor x_dims = 1:x_sz
                if(viz)
                [e j k_dims x_dims]
                end
                x = (2^j)*X(x_dims,:)' - K(k_dims,:)';
                [phis,psis] = arrayfun(@(t) MeyerWavelet(t),x);
                phis = (2^j)*phis';
                psis = (2^j)*psis';
                phi_e = (1-E(e,:)).*phis;
                psi_e = E(e,:).*psis;
                bigPsi_temp = phi_e+psi_e;
                bigPsi_ejk_x = prod(bigPsi_temp);
                % check if vectorization preserves locations in this %
                g_x = dist_array_xy(x_dims);
                summation_C = summation_C + (dx^dim)*g_x*bigPsi_ejk_x; %check use of dim and dx here
            end
            C_ejk(e,j_ix,k_dims) = summation_C;
        end
    end
    
end



