function disturbance = sim_disturbance(y1,m1,d1,y2,m2,d2)

load ClimDat1971_2000.mat

chan=[1 2 3];

n1=datenum(y1,m1,d1);
n2=datenum(y2,m2,d2);
n3=datenum(y1,1,1);

i=find(ClimDat(:,1)>=n1 &  ClimDat(:,1)<=n2); %#ok<NODEF>

t=24*3600*(ClimDat(i,1)-ClimDat(i(1),1));
Data=ClimDat(i,1+chan);

% ClDa=[t Data];
% set_param(gcs,'StopTime',num2str(max(t)));
%set_param([gcs '/StartSec'],'Constant value',num2str(min(t)));

StartSec=24*3600*(ClimDat(i(1),1)-n3);
time = t+StartSec;
sol_diff_hoz = Data(:,1);
sol_dir_norm = Data(:,3);
Qsolar = zeros(length(sol_diff_hoz),1);

for idx = 1:length(sol_diff_hoz)
    Az = 0;
    Inc = 90;
    WindOpp = 1;
    WindZTA = 0.7;
    Qwin1 = irradwindow(time(idx), sol_diff_hoz(idx), sol_dir_norm(idx), Az, Inc, WindOpp, WindZTA);
    
    Az = 90;
    Inc = 90;
    WindOpp = 1;
    WindZTA = 0.8;
    Qwin2 = irradwindow(time(idx), sol_diff_hoz(idx), sol_dir_norm(idx), Az, Inc, WindOpp, WindZTA);
    
    Az = 270;
    Inc = 90;
    WindOpp = 1;
    WindZTA = 0.7;
    Qwin3 = irradwindow(time(idx), sol_diff_hoz(idx), sol_dir_norm(idx), Az, Inc, WindOpp, WindZTA);
    
    Az = 0;
    Inc = 0;
    WindOpp = 0.5;
    WindZTA = 0.7;
    Qwin4 = irradwindow(time(idx), sol_diff_hoz(idx), sol_dir_norm(idx), Az, Inc, WindOpp, WindZTA);
    
    Qsolar(idx) = Qwin1+Qwin2+Qwin3+Qwin4;
    
end

localTime = round(rem((t/3600),24));
Qinternal = 500*(localTime>8 & localTime <18);

disturbance.t = round(t);
disturbance.d1 = Data(:,2);
disturbance.d2 = Qinternal;
disturbance.d3 = Qsolar;

