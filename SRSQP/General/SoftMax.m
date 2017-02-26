function [x_max,C] = SoftMax(vec_x,C)
%x_max = SoftMax(vec_x,C)
% also decimates C if C is large that matlab says log(sum(exp(C*vec_x))) is
% Inf
if(nargin<2)
   C = 100; 
end

%%
x_max = (1/C)*log(sum(exp(C*vec_x)));
while(isinf(x_max))
    'decimating C'
   %C=max(C/10,1);
   C = C/10;
   x_max = (1/C)*log(sum(exp(C*vec_x)));
   C;
end