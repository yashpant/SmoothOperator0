function fhat_x = alt_getWavApprox_vector(x,C_00k,D_ejk,k_min,k_max,j_min,j_max,E_dash)

%%
E = E_dash(2:end,:);
dim = log2(size(E,1)+1);
K = permn(k_min:k_max,dim);
k_sz = size(K,1);
summation = 0;
%%
parfor k_dims = 1:k_sz
    k_dims;
    arg_x = (2^0)*x - K(k_dims,:)';
    [phis,~] = arrayfun(@(t) MeyerWavelet(t),arg_x); %for all elements in vector x
    big_Phi = prod(phis);
    summation = summation + C_00k(1,1,k_dims)*big_Phi;
end

%%
for e = 1:size(E,1)
    j_ix = 0;
    for j=j_min:1:j_max
        j_ix = j_ix + 1;
        parfor k_dims = 1:k_sz
            arg_x = (2^j)*x - K(k_dims,:)';
            [phis,psis] = arrayfun(@(t) MeyerWavelet(t),arg_x); %for all elements in vector x
            phi_e = (1-E(e,:)).*phis';
            psi_e = E(e,:).*psis';
            bigPsi_temp = phi_e+psi_e;
            bigPsi_ejk_x = prod(bigPsi_temp);
            summation = summation + D_ejk(e,j_ix,k_dims)*bigPsi_ejk_x;
            
        end
    end
end
fhat_x = summation;