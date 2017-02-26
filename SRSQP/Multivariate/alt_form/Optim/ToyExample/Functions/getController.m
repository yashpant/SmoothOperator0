function u_opt = getController(u_0,U_intersect_set,options,optParams)

 u_opt = fmincon(@(u)main_objfun2_u_toy_using_mex(u,optParams),u_0,U_intersect_set.A,U_intersect_set.b,[],[],[],[],[],options);
