function [r_P,r_P_der] = robustness_always_SafeDist(traj,d_min,need_derivative)
% robustness for always_I norm(p1-p2)^2>=dmin^2
%trajectory Sx (n*N) (n = sys dim, N = steps)
%%
% if(nargin==2)
%    need_derivative = 0;
% end
r_P_der = zeros(6*2*size(traj,2),1);
der_t = zeros(2*6,size(traj,2)); %change this hardwiring
r_P_der_temp = zeros(2*6*size(traj,2),size(traj,2));
%%
dists_sqd = zeros(size(traj,2),1);

for t = 1:size(traj,2) %for all time steps
    p1 = traj(1:3,t);
    p2 = traj(4:6,t);
    dists_sqd(t) = norm(p1-p2)^2;
    if(need_derivative)
        der_t(:,t) = [zeros(3,1);2*(p1-p2);zeros(3,1);-2*(p1-p2)];
        start = 1+(t-1)*6*2;
        finish = t*6*2;
        r_P_der_temp(start:finish,t) = der_t(:,t);
    end
end
Sd = dists_sqd-d_min^2; %robustness = min_I ( norm(p1-p2)^2 - d_min^2)
[r_P,C] = SoftMin(Sd);
%%
if(need_derivative)
    num = 0;
    for i = 1:size(traj,2)
        num = num + exp(-C*(dists_sqd(i)-d_min^2))*r_P_der_temp(:,i);
    end
    den = sum(exp(-C*Sd));
    r_P_der = num/den;
end

