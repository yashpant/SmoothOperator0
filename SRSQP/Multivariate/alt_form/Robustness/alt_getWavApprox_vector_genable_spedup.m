function fhat_x = alt_getWavApprox_vector_genable_parallel(x,C_00k,D_ejk,k_min,k_max,j_min,j_max,E_dash)

%%
E = E_dash(2:end,:);
dim = log2(size(E,1)+1);

K = zeros((numel(k_min:k_max))^dim,dim);

for i = 1:(numel(k_min:k_max))^dim
    K(i,:) = permn2(k_min:k_max,dim,i);
end

k_sz = size(K,1);
summation = 0;

phis = zeros(1,dim);
psis = zeros(1,dim);

%% get C


for k_dims = 1:k_sz
    arg_x = (2^0)*x - K(k_dims,:)';
    for dimen = 1:dim
        [phis(dimen)] = MeyerWavelet_phi(arg_x(dimen));
    end
    big_Phi = prod(phis);
    summation = summation + C_00k(1,1,k_dims)*big_Phi;
end


%% get D
for e = 1:size(E,1)
    j_ix = 0;
    for j=j_min:1:j_max
        j_ix = j_ix + 1;
        for k_dims = 1:k_sz
            arg_x = (2^j)*x - K(k_dims,:)';
            
            
            for dimen = 1:dim
                [phis(dimen),psis(dimen)] = MeyerWavelet(arg_x(dimen));
            end
            phi_e = (1-E(e,:)).*phis;
            psi_e = E(e,:).*psis;
            
            bigPsi_temp = phi_e+psi_e;
            bigPsi_ejk_x = prod(bigPsi_temp);
            summation = summation + D_ejk(e,j_ix,k_dims)*bigPsi_ejk_x;
            
        end
    end
end
fhat_x = summation;