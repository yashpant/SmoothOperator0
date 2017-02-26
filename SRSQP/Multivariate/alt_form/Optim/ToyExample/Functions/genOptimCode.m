function [status] = genOptimCode(x_0,optParams)


%%
    cfg=coder.config('mex');
        
    %for get coeffs (grid_x,dist_array_xy,dx,j_min,j_max,k_min,k_max,E_dash,0);
    arg_ins = {coder.typeof(x_0),coder.typeof(optParams)};
   %coder.varsize('grid_x',[1 10^5],[0 1]);
    codegen -config cfg objfun2_toy -report -args arg_ins
    codegen -config cfg confun2_toy -report -args arg_ins
    
    status = 1;
end
        
