function bl_r = getBldgRobustness_exact(room_temp,optParams,I1)
%bldg get shizzle
wavparams = optParams.wavparams;
%%

%I2 if needed
room_temp_occu = room_temp(I1);

signed_dists1 = zeros(numel(I1),1);

for i = 1:numel(I1)
    
    x = room_temp_occu(i); %scale here
    signed_dists1(i) = SignedDist(x,optParams.P1_comfort.A,optParams.P1_comfort.b);
         
    
end

bl_r = min(signed_dists1);