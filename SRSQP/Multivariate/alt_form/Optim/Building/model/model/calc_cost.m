function cost = calc_cost(u)

ethaHeat=0.9;
ethaCool=3;
costHeat=0.25;
costCool=0.2;
UnitHeat=1;
UnitCool=1;

cost = buildingv5efun([ethaHeat, ethaCool, costHeat, costCool, UnitHeat, UnitCool, u]);