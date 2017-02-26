function [C,D] = getCoefficientsScalar(grid_x,dist_array_x,dx,j_min,j_max,k_min,k_max,visualization)
% [C,D] = getCoefficientsScalar(grid_x,dist_array_x,dx,j_min,j_max,k_min,k_max,visualization)

%% get coeffs c's and d's, as in Wavelets, Phillip K. Poon
%visualization  = 1;

j_ix = 0;
k_ix = 0;

C = zeros(numel(k_min:1:k_max),1);
D = zeros(numel(j_min:1:j_max), numel(k_min:1:k_max));

for j=j_min:1:j_max
    j_ix = j_ix + 1;
    k_ix = 0;
    for k=k_min:1:k_max
        k_ix = k_ix+1;
        summation_c = 0;
        summation_d = 0;
        for ix = 1:numel(grid_x)
            x = grid_x(ix);
            g = dist_array_x(ix); %function
            
            [phi_jk,psi_jk] = MeyerWavelet(2^j*x-k); %scale and wav
            phi_jk = (2^(j/2))*phi_jk;
            psi_jk = (2^(j/2))*psi_jk;
            if(j==0) %coeffs for phi (note, j=0 is base scale);
                summation_c = summation_c + dx*g*phi_jk;
            end
            %coeffs for psi
            summation_d = summation_d + dx*g*psi_jk;
        end
        if(j==0)
            C(k_ix) = summation_c;
        end
        D(j_ix,k_ix) = summation_d;
    end
end

%%
if(visualization)
    figure;
    plot(k_min:1:k_max,C);grid on;xlabel('k');ylabel('c(0,k)');
    figure;
    if(j_min==j_max)
        plot(k_min:1:k_max,D);
        xlabel('k');
    else
        mesh(k_min:1:k_max,j_min:1:j_max,D);grid on;
        xlabel('k');
        ylabel('j');
    end
    
end