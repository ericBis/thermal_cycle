function [DAT_WATER DAT_AIR MASSFLOW] = CoolingTower(P_w,options)
% COOLINGTOWER is a cooling tower 0D modelisation
% COOLINGTOWER(P_w,options) compute the thermodynamics states for a Cooling
% tower based on several inputs (given in OPTION) and based on a given 
% water power to dissipate P_w.
% It returns the main results. 
%
% INPUTS :
% P_W = Heat power output at the condenser [kW]
% OPTIONS is a structure containing :
%   -options.Tcond  [°C]: Temperature in the condenser
%   -options.Tpinch [°C]: Minimum tempearture pinch between Tw_out and the
%                         condenser temperature.
%   -options.Tw_out [°C]: Cooling water temperature at the condenser outlet
%   -options.Tw_in  [°C]: Cooling water temperature at the condenser inlet
%   -options.Triver [°C]: River temperature 
%   -options.Ta_in  [°C]: Atmospheric air temperature 
%   -options.Ta_out [°C]: Air outlet temperature of cooling tower 
%   -options.Phi_atm [-]: Relative humidity of atmospheric air
%   -options.Phi_out [-]: Maximum relative humidity of air at the cooling 
%                         tower outlet.
%
% OUTPUT :
% MassFlow [kg/s]: Vector containing the different massflow :
%   -massflow(1) : water massflow at the condenser
%   -massflow(2) : additionnal water massflow = water flow evaporated
%   -massflow(3) : air massflow at the cooling tower   
%
%  dat_water = [T_e1       , T_e2       , T_e3       , T_e4;  %[°C]
%               h_e1       , h_e2       , h_e3       , h_e4;  %[kJ/kg]
%               m_e1       , m_e2       , m_e3       , m_e4]; %[kg/s]
% 
%  dat_air   = [Ta_in       , Ta_out  ;  %[°C]
%               ha_in       , ha_out  ;  %[kJ/kg]
%               xa_in       , xa_out  ;  %[kg_water/kg_dry_air]
%               Phia_in     , Phia_out]; %[-] relative humidity
%  
%
% ADDITIONNAL INFORMATIONS
% Water points : 
%       1 : water outlet of cooling tower
%       2 : water just before condenser
%       3 : water just after  condenser
%       4 : water from the river (coming between 1 & 2)
%
% Air points :
%       a_in : air at the cooling tower inlet
%       a_out : air at the cooling tower outlet
%

%% YOUR WORK
if nargin<2
    options=struct();
    if nargin<1
        P_w=200e3;%200MW_heat
    end
end

if isfield(options,'Tpinch')
    Tpinch = options.Tpinch;
else
    Tpinch = 4; %[K]
end

DAT_AIR=zeros(4,2); % Prealocations 
DAT_WATER=zeros(3,4);
MASSFLOW=zeros(3,1);

Tcond=options.Tcond;
Tpinch=options.Tpinch;
Tw_out=options.Tw_out;
Tw_in=options.Tw_in;
Triver=options.Triver;
Ta_in=options.Ta_in;
Ta_out=options.Ta_out;
Phi_atm=options.Phi_atm;
Phi_out=options.Phi_out;

te2=Ta_in;
te3=Tw_out;
te4=Triver;

DAT_AIR(1,1)=Ta_in;
DAT_AIR(4,1)=Phi_out;
[~, xa_in, ~, ha_in, ~, ~, ~] = Psychrometrics('Tdb',Ta_in,'phi',Phi_atm*100);
DAT_AIR(2,1)=ha_in;
DAT(3,1)=xa_in;

he2=XSteam('h_Tx',te2,0);
he3=XSteam('h_Tx',te3,0);
he4=XSteam('h_Tx',te4,0);

eq = @(x) Eq_CL(x,Phi_out,ha_in,he2*1e3,he3*1e3,he4*1e3,xa_in,P_w*1e3);

x0= [2*ha_in ; 10*xa_in ; 20; 38 ; 40 ; 1 ; 0.6*he3]; 

x = fsolve(eq,x0);




end