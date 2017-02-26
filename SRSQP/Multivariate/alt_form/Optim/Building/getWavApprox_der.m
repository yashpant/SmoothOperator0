function fhat_dash_x = getWavApprox_der(x,C,D,k_min,k_max,j_min,j_max)
%% fhat_x = getWavApprox(x,C,D,k_min,k_max,j_min,j_max)

sum_1 = 0;
sum_2 = 0;
k_ix = 0;
j_ix = 0;
KJ = [];
for k = k_min:1:k_max
    k_ix = k_ix + 1;
    j_ix = 0;
    [phi_xminusk,~] = MeyerWavelet_der(x-k);
    sum_1 = sum_1 + C(k_ix)*phi_xminusk;
    for j = j_min:1:j_max
        j_ix = j_ix + 1;
        KJ = [KJ;k_ix j_ix]; %sanity check
        [~,psi_2jmk] = MeyerWavelet_der((2^j)*x-k);
       sum_2 = sum_2 + D(j_ix,k_ix)*psi_2jmk;
    end
end

fhat_dash_x = sum_1 + sum_2;



