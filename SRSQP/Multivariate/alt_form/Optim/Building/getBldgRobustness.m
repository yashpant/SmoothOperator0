function [bl_r,bl_r_der] = getBldgRobustness(room_temp,optParams,I1,der)
%bldg get shizzle in interval I1
bl_r_der = [];
wavparams = optParams.wavparams;
%%

%I2 if needed
room_temp_occu = room_temp(I1);

signed_dists1 = zeros(numel(I1),1);
signed_dists_der = zeros(numel(I1),1);

parfor i = 1:numel(I1)
    x = room_temp_occu(i)/10; %scale here
    signed_dists1(i) = getWavApprox(x,wavparams.C,wavparams.D, ...
        wavparams.k_min,wavparams.k_max,wavparams.j_min,wavparams.j_max);
end

[bl_r,C] = SoftMin(10*signed_dists1);



if(der)
    bl_r_der = zeros(numel(room_temp),1);
    dSmin_by_du = exp(-C*10*signed_dists1)./sum(exp(-C*10*signed_dists1));
    parfor i = 1:numel(I1) %all time instants in consideratiojn
        x = room_temp_occu(i)/10; %scale here
        signed_dists_der(i) = getWavApprox_der(x,wavparams.C,wavparams.D, ...
            wavparams.k_min,wavparams.k_max,wavparams.j_min,wavparams.j_max);
    end
    bl_r_der(I1) = dSmin_by_du.*signed_dists_der;
end

