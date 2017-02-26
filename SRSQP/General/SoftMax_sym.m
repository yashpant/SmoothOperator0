function x_max = SoftMax_sym(vec_x,C)
%x_max = SoftMax(vec_x,C)
% also decimates C if C is large that matlab says log(sum(exp(C*vec_x))) is
% Inf
if(nargin<2)
   C = 100; 
end
if(isa(vec_x,'sym'))
   syms C; 
end

%%
x_max = (1/C)*log(sum(exp(C*vec_x)));
if(~isa(vec_x,'sym'))
while(isinf(x_max))
    'decimating C'
   C=max(C/10,1);
   x_max = (1/C)*log(sum(exp(C*vec_x)))
   C;
end
end