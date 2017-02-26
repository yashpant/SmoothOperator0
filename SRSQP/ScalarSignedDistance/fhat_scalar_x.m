function fhat = fhat_scalar_x(x,theta_jk,dist_array_x,dx,J,K)
%fhat = fhat_scalar_x(x,theta_jk,dist_array_x,dx,J,K)

%% evaluate c_0 = integral_-inf^inf f(t)dt
c_0 = dx*sum(dist_array_x(:));

%% evaluate f_hat
fhat = 0;%c_0;

for j_ix = 1:numel(J)
    for k_ix = 1:numel(K)
        j = J(j_ix);
        k = K(k_ix);
        arg_x = ((2^j)*x)-k;
        if(0)
        psi_jk = (2^(1*j/2))*RickerWavelet(arg_x);
        else
        [phi,psi] = MeyerWavelet(arg_x);
        psi_jk = (2^(1*j/2))*psi;
        end
        fhat = fhat + theta_jk(j_ix,k_ix)*psi_jk;
    end
end