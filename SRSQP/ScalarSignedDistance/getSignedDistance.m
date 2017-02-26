function [signed_distance,close_pt] = getSignedDistance(x,P)
% compute signed distance between polytope P and point x, requires 
% MPT and CVX

% input:
% P: MPT Polyhedron
% x: point to compute signed distance for
% output:
% signed_distance
% close_pt: one minimizer of that dist in P (or its edges)


dim = size(P.A,2);
numConstrs = size(P.A,1);
optvals = zeros(numConstrs,1);
optpts = zeros(dim,numConstrs);

if(isInside(P,x))
    for i = 1:numConstrs
        cvx_clear
        cvx_begin quiet
        variable y(dim,1)
        minimise norm(y-x)
        subject to
        P.A(i,:)*y == P.b(i);
        cvx_end
        optvals(i) = cvx_optval;
        optpts(:,i) = cvx_optpnt.y;
    end
    [signed_distance,ix] = min(optvals);
    close_pt = optpts(:,ix);
    
else %if x not in P
    cvx_clear
    cvx_begin quiet
    variable y(dim,1)
    minimise norm(y-x)
    subject to
    P.A*y <= P.b;
    cvx_end
    signed_distance = -cvx_optval; %negative outside
    close_pt = cvx_optpnt.y;
end