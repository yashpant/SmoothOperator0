% example 2
%(always not P1) AND/OR (min/max) (always not P2)
P1 = Polyhedron('A',1,'b',-2);
P2 = Polyhedron('A',-1,'b',-2);
[Params_P1] = WavSignedDistScalar(P1,xmin,xmax,dx,0);
[Params_P2] = WavSignedDistScalar(P2,xmin,xmax,dx,0);
%pause;
close all;clc;
%plot(P1);hold all;plot(P2);

%% 
Ntrials = 100;

for i = 1:Ntrials

trajec(:,i) = 0.5*randn(10,1);
traj = trajec(:,i);

r_exact_P1 = -getRobustnessP(traj,P1,Params_P1,1);
r_approx_P1 = -getRobustnessP(traj,P1,Params_P1,0);

r_exact_P2 = -getRobustnessP(traj,P2,Params_P2,1);
r_approx_P2 = -getRobustnessP(traj,P2,Params_P2,0);

r_phi_exact = min(r_exact_P1,r_exact_P2);
r_phi_approx = SoftMin([r_approx_P1,r_approx_P2]);
err(i) = r_phi_exact-r_phi_approx;
err_rel(i) = (r_phi_exact-r_phi_approx)/abs(r_phi_exact);
end
%[mean(abs(err)) std(abs(err))]
[mean(err_rel) std(err_rel)]
figure;
hist(err_rel);grid on;