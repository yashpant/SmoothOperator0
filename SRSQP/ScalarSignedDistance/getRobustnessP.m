function r_P = getRobustnessP(traj,P,wavparams,exact)
%r_P = getRobustnessP(traj,P,wavparams,exact)
%get robustness of a scalar system's , w.r.t a Polynomial P as the
%predicate set

if(nargin~=4)
    exact = 0;
end
if(nargin<3)
   'come on man ...'
   return
end

% if exact
if(exact==1)
    %get signed distances
    if(exist('Staliro')>0)
        signed_dists = arrayfun(@(x) SignedDist(x,P.A,P.b), traj);
    else
        signed_dists = arrayfun(@(x) getSignedDistance(x,P), traj);
    end
    
    %minimum over instances
    dmin = min(signed_dists);
    
    %else if approx
else
    signed_dists = arrayfun(@(x) getWavApprox(x,wavparams.C,wavparams.D,...
        wavparams.k_min,wavparams.k_max,wavparams.j_min,wavparams.j_max),...
        traj);
    dmin = SoftMin(signed_dists);
end

r_P = dmin;
