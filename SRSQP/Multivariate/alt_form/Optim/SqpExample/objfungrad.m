function [f,gradf] = objfungrad(x)
f = max(x(1),x(2));
%f = exp(x(1))*(4*x(1)^2+2*x(2)^2+4*x(1)*x(2)+2*x(2)+1);
% Gradient of the objective function:
if nargout  > 1
    gradf = [ f + exp(x(1)) * (8*x(1) + 4*x(2)), 
    exp(x(1))*(4*x(1)+4*x(2)+2)];
end