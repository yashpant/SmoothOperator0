function [r_P,Sd] = getRobustnessP_vector(traj,P,wavparams,exact)
%trajectory Sx (n*N) (n = sys dim, N = steps)
%%
S = num2cell(traj,1)';
 
if(~exact)

signed_dists_inex = cellfun(@(x) getWavApprox_vector(x,wavparams.C_ejk,...
        wavparams.k_min,wavparams.k_max,wavparams.j_min,wavparams.j_max, ...
        wavparams.E), S, 'UniformOutput', false);
C = signed_dists_inex;
else
    if(exist('Staliro')>0)
    signed_dists = cellfun(@(x) SignedDist(x,P.A,P.b), S, 'UniformOutput', false);
    
    else
    signed_dists = cellfun(@(x) getSignedDistance(x,P), S, 'UniformOutput', false);
    end
    C = signed_dists;
end

S1 = sprintf('%s*', C{:});
Sd = sscanf(S1, '%f*');

r_P = min(Sd);