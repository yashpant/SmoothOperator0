function f = objfun_feas(x,optParams)

f = norm(x(1:optParams.dim*optParams.len),2);
