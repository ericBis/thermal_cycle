Pw=187297; % [kW]
options= struct;
options.Tcond= 35; % [°C]: Temperature in the condenser
options.Tpinch = 0;% [°C]:Minimum tempearture pinch between Tw_out and the
%                         condenser temperature.
options.Tw_out = 22.5;%[°C] Cooling water temperature at the condenser outlet
options.Tw_in = 14.2; %[°C] Cooling water temperature at the condenser inlet
options.Triver= 15;% [°C] River temperature
options.Ta_in = 5;% [°C] Atmospheric air temperature 
options.Ta_out= 16;% [°C] Air outlet temperature of cooling tower
options.Phi_atm= 0.45;% [-] Relative humidity of atmospheric air
options.Phi_out= 1;% [-] Maximum relative humidity of air at the cooling 
%                         tower outlet.

[DAT_WATER DAT_AIR MASSFLOW] = CoolingTower(Pw,options);
