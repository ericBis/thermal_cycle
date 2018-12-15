function [DAT_WATER, DAT_AIR, MASSFLOW] = CoolingTower(P_w,options)
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

if ~isfield(options,'Tw_out')
    options.Tw_out=60; % [°C]
end

if ~isfield(options,'Tw_in')
    options.Tw_in=30; % [°C]
end

if ~isfield(options,'Triver')
    options.Triver=15; % [°C]
end

if ~isfield(options,'Ta_in')
    options.Ta_in=15; % [°C]
end

if ~isfield(options,'Ta_out')
    options.Ta_out=15; % [°C]
end

if ~isfield(options,'Phi_atm')
    options.Phi_atm=0.6;
end

if ~isfield(options,'Phi_out')
    options.Phi_out=0.99; 
end

%Tcond=options.Tcond;
Tw_out=options.Tw_out;
Tw_in=options.Tw_in;
Triver=options.Triver;
Ta_in=options.Ta_in;
Ta_out=options.Ta_out;
Phi_atm=options.Phi_atm;
Phi_out=options.Phi_out;
%p_atm=1.01325; % [bar] pressure of the atmosphere

% Air caractéristics
[~, xa_in, ~, ha_in, ~, ~, ~] = Psychrometrics('Tdb',Ta_in,'phi',Phi_atm*100);
[~, xa_out, ~, ha_out, ~, ~, ~] = Psychrometrics('Tdb',Ta_out,'phi',Phi_out*100);

%water caracteristics
%te1=Triver+Tpinch;
te2=Tw_in;
te3=Tw_out;
te4=Triver;
%he1=XSteam('h_pT',p_atm,te1);
he2=XSteam('hL_T',te2);
he3=XSteam('hL_T',te3);
he4=XSteam('hL_T',te4);

%massflows
m_as=P_w*1e3/(ha_out-ha_in);
m_evap=m_as*(xa_out-xa_in);
me2=P_w/(he3-he2);
me3=me2;
me4=m_evap;
me1=me3-m_evap;

%water state 1 caracteristics
he1=(me2*he2-me4*he4)/me1;
te1=XSteam('T_ph',1,he1);

% Outputs
DAT_AIR = [Ta_in Ta_out ; ha_in ha_out ; xa_in xa_out ; Phi_atm Phi_out];
DAT_WATER = [te1 te2 te3 te4 ; he1 he2 he3 he4 ; me1 me2 me3 me4];
MASSFLOW = [me2 ; m_evap ; m_as];

end