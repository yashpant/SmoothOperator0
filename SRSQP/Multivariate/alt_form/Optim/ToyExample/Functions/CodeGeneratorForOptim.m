% wavparams
wavparams = optParams.Params_P_unsafe;
traj = reshape(x_0(1:optParams.dim*optParams.len),optParams.dim,optParams.len);
P = [];
%% codegen for
%alt_getRobustnessP_vector_genable_parallel(traj,P,wavparams,exact)

if(1)
cfg=coder.config('mex');
arg_ins = {coder.typeof(traj),coder.typeof(P),coder.typeof(wavparams),coder.typeof(0)};
%coder.varsize('grid_x',[1 10^5],[0 1]);
codegen -config cfg alt_getRobustnessP_vector_genable_parallel -report -args arg_ins
end
%%
if(0)
%parallel was faster
codegen -config cfg alt_getRobustnessP_vector_genable -report -args arg_ins
end
%%
%alt_getRobustnessP_and_der_vector_genable(traj,wavparams)
if(1)
arg_ins = {coder.typeof(traj),coder.typeof(wavparams),coder.typeof(0)};
codegen -config cfg alt_getRobustnessP_and_der_vector_genable -report -args arg_ins
end

%% parallel
if(1)
arg_ins = {coder.typeof(traj),coder.typeof(wavparams),coder.typeof(0)};
codegen -config cfg alt_getRobustnessP_and_der_vector_genable_parallel -report -args arg_ins
end
%% for eventually
%robustness_eventually_P_genable(traj,wavparams,need_derivative)
if(0)
arg_ins = {coder.typeof(traj),coder.typeof(wavparams),coder.typeof(0)};
codegen -config cfg robustness_eventually_P_genable -report -args arg_ins
end
%% parallel
if(1)
arg_ins = {coder.typeof(traj),coder.typeof(wavparams),coder.typeof(0)};
codegen -config cfg robustness_eventually_P_genable_parallel -report -args arg_ins
end

%% objfun2_toy_genable(x,optParams)
if(0)
arg_ins = {coder.typeof(x_0),coder.typeof(optParams),coder.typeof(1)};
codegen -config cfg objfun2_toy_genable -report -args arg_ins
end

%%
% cfg=coder.config('mex');
%   arg_ins = {coder.typeof(u_0),coder.typeof(U_intersect_set),coder.typeof(options),coder.typeof(optParams)}
% codegen -config cfg getController -report -args arg_ins
