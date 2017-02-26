%%
%wavparams = SmoothOpt.preds.WavParams;
wavparams(1) = optParams.Params_P_unsafe;
wavparams(2) = optParams.Params_P_term;
%%
x = [0;0]
%x = [1;1]
alt_getWavApprox_vector_genable(x,wavparams(1).C_00k,wavparams(1).D_ejk,wavparams(1).k_min, ...
    wavparams(1).k_max,wavparams(1).j_min,wavparams(1).j_max,wavparams(1).E_dash)

%%
x = [0;0];
x = [2.25;2.25];
x = traj_x(:,end);

alt_getWavApprox_vector_genable(x,wavparams(2).C_00k,wavparams(2).D_ejk,wavparams(2).k_min, ...
    wavparams(2).k_max,wavparams(2).j_min,wavparams(2).j_max,wavparams(2).E_dash)