function [r_P,Sd] = alt_getRobustnessP_vector_in(traj,P,wavparams,exact)
%trajectory Sx (n*N) (n = sys dim, N = steps)
%%
S = num2cell(traj,1)';
 
if(~exact)

signed_dists_inex = cellfun(@(x) alt_getWavApprox_vector(x,wavparams.C_00k, wavparams.D_ejk,...
        wavparams.k_min,wavparams.k_max,wavparams.j_min,wavparams.j_max, ...
        wavparams.E_dash), S, 'UniformOutput', false);

C = signed_dists_inex;
format long
Sd = +cell2mat(C); %check + - based on formula
%keyboard
%S1 = sprintf('%f*', C{:});
%Sd = sscanf(S1, '%f*');
r_P = SoftMin(Sd);
else
    if(exist('Staliro')>0)
    signed_dists = cellfun(@(x) SignedDist(x,P.A,P.b), S, 'UniformOutput', false);
    
    else
    signed_dists = cellfun(@(x) getSignedDistance(x,P), S, 'UniformOutput', false);
    end
   C = signed_dists;
Sd = +cell2mat(C);
  
%S1 = sprintf('%f*', C{:});
%Sd = sscanf(S1, '%f*');

r_P = min(Sd);
end


