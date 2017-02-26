%BUILDV4EFUN 	Building simulatie Analyse
%				
%JvS 2001,sept


function cost=buildingv5efun(x)

ethaHeat=x(1);
ethaCool=x(2);
costHeat=x(3);
costCool=x(4);
UnitHeat=x(5);
UnitCool=x(6);
Q=x(7);

y=zeros(4,1);

H=32e6;  %J/m3
kWh=3.6e6; %J/kWh


if Q>0
y(1)=Q/(ethaHeat*H);  
else
y(2)=-Q/(ethaCool*H);    
end


if UnitHeat==1 %[Eur/m3]
   y(3)=y(1)*costHeat;
else
   y(3)=(H/kWh)*y(1)*costHeat;
end

if UnitCool==1 %[Eur/m3]
   y(4)=y(2)*costCool;
else
	y(4)=(H/kWh)*y(2)*costCool;   
end

cost = y(3)+y(4);



