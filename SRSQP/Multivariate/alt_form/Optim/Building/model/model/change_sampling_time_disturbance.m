function d = change_sampling_time_disturbance(dist, time)

% sample disturbance on the same time scale using interpolation or keep it
% constant for 1hr duration?
ifInerpolateDist = 0;
d = zeros(3,length(time));
if ifInerpolateDist
    d(1,:) = interpc(dist.t, dist.d1, time); %#ok<*UNRCH>
    d(2,:) = interpc(dist.t, dist.d2, time);
    d(3,:) = interpc(dist.t, dist.d3, time);
else
    idx = 1;
    for idt = 1:length(time)
        while idx<=length(dist.t) && dist.t(idx)<=time(idt)
            idx = idx+1;
        end
        d(1,idt) = dist.d1(idx-1);
        d(2,idt) = dist.d2(idx-1);
        d(3,idt) = dist.d3(idx-1);
    end
end