%IRRADWINDOW  Calculates the irradiation through window
%
% SYNTAX: E=irrad(Dh,En,iday,LST,gamma,beta)
%
% OUTPUT: vector E  
% E(1)= diffuse solar irradiation on an inclined surface
% E(2)= direct solar irradiation on an inclined surface
% E(3)= total solar irradiation on an inclined surface
% E(4)= total solar irradiation on a horizontal surface
%
% INPUT:
% (scalar) Dh    = diffuse horizontal irradiation     [W/m2]
% (scalar) En    = direct normal irradiation          [W/m2]
% (scalar) t     =  time seconds after midnight 1 january
% (scalar) gamma = azimuth angle of the surface,   
%             east:gamma = -90,     west:gamma = 90
%             south:gamma = 0,    north:gamma = 180
% (scalar) beta  = inclination angle of the surface,
%             horizontal: beta=0, vertical: beta=90
%
% default geographical position: De Bilt
% default ground reflectivity (albedo): 0.2
%
% EXAMPLE: E=irrad(800,200,201,12,0,45)
% ANSWER: E=1.0e+003 *
%    0.8569    0.1907    1.0759    0.9684
%
% REF: Perez (zie Solar Energy volume 39 no. 3)
%
% JvS feb 2002, aanpassing voor inlezen standaard klimaatfiles

function EW=irradwindow(tclim,Dh,En,gamma,beta,WA,WZTA);

%(scalar) iday  = day of the year                   (1-365)    
% (scalar) LST   = Local Standard time (0 - 24)       [hour]

iday=1+floor(tclim/(24*3600));
LST=floor(rem( (tclim/3600),24));


% L   = Latitude [graden] 
L=52.1;
% LON = Local Longitude [graden] oost is positief
LON=5.1;
% LSM = Local Standard time Meridian [graden] oost is positief
LSM=15;
% gref = albedo
gref=0.2;


r=pi/180;
L=L*r;
beta=beta*r;
theta=2*pi*(iday-1)/365.25;
el=4.901+0.033*sin(-0.031+theta)+theta;
% declination
delta=asin(sin(23.442*r)*sin(el));
q1=tan(4.901+theta);
q2=cos(23.442*r)*tan(el);
% equation of time
ET=(atan((q1-q2)/(q1*q2+1)))*4/r;
AST=LST+ET/60-(4/60)*(LSM-LON);
h=(AST-12)*15*r;

% hai=sin(solar altitude)
hai=cos(L)*cos(delta)*cos(h)+sin(L)*sin(delta);

E(1)=0; E(2)=0; E(3)=0; E(4)=0;
if hai>0,


% salt=solar altitude
salt=asin(hai);
phi=acos((hai*sin(L)-sin(delta))/(cos(salt)*cos(L)))*sign(h);
gam=phi-gamma*r;

% cai=cos(teta)
cai=cos(salt)*cos(abs(gam))*sin(beta)+hai*cos(beta);
% teta = incident angle on the tilted surface
teta=acos(cai);
% salts=solar altitude for an inclined surface
salts=pi/2-teta;

% Perez (zie Solar Energy volume 39 no. 3)
% berekening van de diffuse straling op een schuin vlak
% Approximatin of A and C, the solid angles occupied by the circumsolar region,
% weighed by its average incidence on the slope and horizontal respectively.
% In the expression of diffuse on inclined surface the quotient of A/C is
% reduced to XIC/XIH. A=2*(1-cos(beta))*xic, C=2*(1-cos(beta))*xih
% gecontroleerd  okt 1996 martin de wit

% alpha= the half-angle circumsolar region
alpha=25*r;

if salts<-alpha,
   xic=0;
elseif salts>alpha,
   xic=cai;
else
   xic=0.5*(1+salts/alpha)*sin((salts+alpha)/2);
end

if salt>alpha,
   xih=hai;
else
   xih=sin((alpha+salt)/2);
end

epsint=[1.056 1.253 1.586 2.134 3.23 5.98 10.08 999999];
f11acc=[-0.011 -0.038 0.166 0.419 0.710 0.857 0.734 0.421];
f12acc=[0.748 1.115 0.909 0.646 0.025 -0.370 -0.073 -0.661];
f13acc=[-0.080 -0.109 -0.179 -0.262 -0.290 -0.279 -0.228 0.097];
f21acc=[-0.048 -0.023 0.062 0.140 0.243 0.267 0.231 0.119];
f22acc=[0.073 0.106 -0.021 -0.167 -0.511 -0.792 -1.180 -2.125];
f23acc=[-0.024 -0.037 -0.050 -0.042 -0.004 0.076 0.199 0.446];

% determination of zet = solar zenith angle (pi/2 - solar altitude).
zet=pi/2-salt; 

% determination of inteps with eps
inteps=1;
if Dh>0,
   eps=1+En/Dh;
   i=find(epsint>=eps);
   inteps=min(i);
end

% calculation of inverse relative air mass
airmiv=hai;
if salt<10*r,
   airmiv=hai+0.15*(salt/r+3.885)^(-1.253);
end

% calculation of extraterrestrial radiation
Eon=1370*(1+0.033*cos(2*pi*(iday-3)/365));

% delta is "the new sky brightness parameter"
delta=Dh/(airmiv*Eon);

% determination of the "new circumsolar brightness coefficient
% (f1acc) and horizon brightness coefficient (f2acc)"
f1acc=f11acc(inteps)+f12acc(inteps)*delta+f13acc(inteps)*zet;
f2acc=f21acc(inteps)+f22acc(inteps)*delta+f23acc(inteps)*zet;

% determination of the diffuse radiation on an inclined surface
E(1)=Dh*(0.5*(1+cos(beta))*(1-f1acc)+f1acc*xic/xih+f2acc*sin(beta));
if E(1)<0,
   E(1)=0;
end

% horizontal surfaces treated separately
% beta=0 : surface facing up, beta=180(pi) : surface facing down
if beta>-0.0001 & beta<0.0001,
   E(1)=Dh;
end
if beta>(pi-0.0001) & beta<(pi+0.0001),
   E(1)=0;
end

% Isotropic sky
% E(1)=0.5*(1+cos(beta))*Dh;

% direct solar radiation on a surface
E(2)=En*cai;
if E(2)<0.0, 
   E(2)=0;
end

% the ground reflected component: assume isotropic
% ground conditions.
Eg=0.5*gref*(1-cos(beta))*(Dh+En*hai);

% global irradiation
E(4)=Dh+En*hai;

% total irradiation on an inclined surface
E(3)=E(1)+E(2)+Eg;

end

EW=E(3)*WA*WZTA;