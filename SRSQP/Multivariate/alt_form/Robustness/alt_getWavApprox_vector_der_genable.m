function fhat_x = alt_getWavApprox_vector_der_genable(x,ix,C_00k,D_ejk,k_min,k_max,j_min,j_max,E_dash)

%derivative of signed dist approximation w.r.t ix^th element of x

%%
E = E_dash(2:end,:);
dim = log2(size(E,1)+1);

K = zeros((numel(k_min:k_max))^dim,dim);

for i = 1:(numel(k_min:k_max))^dim
    K(i,:) = permn2(k_min:k_max,dim,i);
end

%K = permn(k_min:k_max,dim);
k_sz = size(K,1);
summation = 0;

phis_0 = zeros(1,dim-1);
phis = zeros(1,dim);
psis = zeros(1,dim);
phi_d = zeros(1);
psi_d = zeros(1);

dims = setdiff(1:dim,ix);


%%
for k_dims = 1:k_sz
    arg_x = (2^0)*x - K(k_dims,:)';
    
    if(dim>1) %if not scalar, if scalar, let phis_0 be 0
        for dimen = 1:dim-1
            [phis_0(dimen),~] = MeyerWavelet(arg_x(dims(dimen)));
        end
    end
    %the derivative of Phi_0k(arg_x)
    [phi_d,~] = MeyerWavelet_der(arg_x(ix));
    %[phis,~] = arrayfun(@(t) MeyerWavelet(t),arg_x); %for all elements in vector x
    if(dim>1)
        big_Phi = prod(phis_0)*phi_d;
    else
        big_Phi = phi_d;
    end
    summation = summation + C_00k(1,1,k_dims)*big_Phi;
end
%%

for e = 1:size(E,1)
    j_ix = 0;
    for j=j_min:1:j_max
        j_ix = j_ix + 1;
        for k_dims = 1:k_sz
            arg_x = (2^j)*x - K(k_dims,:)';
            
            if(dim>1)
                for dimen = 1:dim-1
                    [phis(dims(dimen)),psis(dims(dimen))] = MeyerWavelet(arg_x(dims(dimen)));
                end
            end
            [phi_d, psi_d] = MeyerWavelet_der(arg_x(ix));
            phis(ix) = (2^j)*phi_d;
            psis(ix) = (2^j)*psi_d;
            
            %[phis,psis] = arrayfun(@(t) MeyerWavelet(t),arg_x); %for all elements in vector x
            phi_e = (1-E(e,:)).*phis;
            psi_e = E(e,:).*psis;
            
            bigPsi_temp = phi_e+psi_e;
            bigPsi_ejk_x = prod(bigPsi_temp);
            summation = summation + D_ejk(e,j_ix,k_dims)*bigPsi_ejk_x;
            
        end
    end
end
fhat_x = summation;