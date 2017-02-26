function [A,B,Bd,C,sysc] = linear_model(Ts)

model_params;

if N_interne_massa==1;
   
%zwaarbeton   
c_interne_massa=840;
d_interne_massa=0.2;
rho_interne_massa=2500;
   
   
elseif  N_interne_massa==2;  
   
%middelbeton   
c_interne_massa=840;
d_interne_massa=0.1;
rho_interne_massa=1000;
     
else
   
%licht beton   
c_interne_massa=860;
d_interne_massa=0.1;
rho_interne_massa=500;
     
end

%gebouwschil
if N_gebouwschil==1
   
%zwaarbeton
c_gebouwschil=840;
rho_gebouwschil=2500;
d_gebouwschil=0.2;

   
elseif  N_gebouwschil==2  
   
%middelbeton   
c_gebouwschil=840;
rho_gebouwschil=1000;
d_gebouwschil=0.1;

else
   
%lichtbeton   
c_gebouwschil=840;
rho_gebouwschil=500;
d_gebouwschil=0.1;

end


V_interne_massa=A_interne_massa*d_interne_massa;


rho_lucht=1.29;
c_lucht=1005;
alpha_i_gebouwschil=8;
alpha_e_gebouwschil=23;
alpha_i_glas=8;
alpha_e_glas=23;
alpha_interne_massa=8;
deelfactor=0.8;

% formules
fie=(n*V_gebouw)/3600;
m_stroom=fie*rho_lucht;
V_gebouwschil=A_gebouwschil*d_gebouwschil;

% bepaling weerstanden
R1=1/(A_interne_massa*alpha_interne_massa);
R2=1/(A_gebouwschil*alpha_i_gebouwschil);
R3=Rc_gebouwschil_gem/A_gebouwschil;
R4=1/(A_gebouwschil*alpha_e_gebouwschil);
R5=1/(A_glas_werk*alpha_i_glas) + Rc_glas/A_glas_werk + 1/(A_glas_werk*alpha_e_glas);

% bepaling capaciteiten
C1=rho_interne_massa*c_interne_massa*V_interne_massa;
C2=rho_gebouwschil*c_gebouwschil*V_gebouwschil;
C3=C2;
Ci=rho_lucht*c_lucht*V_gebouw;

% bepaling matrix-elementen
a1=-1/(R1*C1);
a2=1/(R1*C1);
a3=-1/(R2*C2)-1/(R3*C2);
a4=1/(R3*C2);
a5=1/(R2*C2);
a6=1/(R3*C3);
a7=-1/(R3*C3)-1/(R4*C3);
a8=1/(R1*Ci);
a9=1/(R2*Ci);
a10=-m_stroom*c_lucht/Ci-1/(R2*Ci)-1/(R5*Ci)-1/(R1*Ci);
b1=deelfactor/C1;
b2=1/(R4*C3);
b3=m_stroom*c_lucht/Ci+1/(R5*Ci);
b4=1/Ci;
b5=1/Ci;
b6=(1-deelfactor)/Ci;

% bepaling matrices
Ac=[a1 0 0 a2 ; 0 a3 a4 a5 ; 0 a6 a7 0 ; a8 a9 0 a10];
Bc=[0 0 0 b1 ; 0 0 0 0 ; b2 0 0 0 ; b3 b4 b5 b6];
Cc=[1 0 0 0 ; 0 1 0 0 ; 0 0 1 0 ; 0 0 0 1];
Dc=[0 0 0 0 ; 0 0 0 0 ; 0 0 0 0 ; 0 0 0 0];

sysc = ss(Ac, Bc, Cc, Dc);
sysd = c2d(sysc, Ts);

A = sysd.a;
B = sysd.b(:,2);
C = sysd.c;
Bd = sysd.b(:,[1,3:4]);
