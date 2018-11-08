%% Simple Rankine-Hirn
P_e=1750; % [kW]
options= struct;
options.nsout=0; %   [-] : Number of feed-heating 
options.reheat=0; %    [-] : Number of reheating
options.T_max=0; %     [�C] : Maximum steam temperature
options.T_cond_out=0; %[�C] : Condenseur cold outlet temperature
options.p3_hp=0; %     [bar] : Maximum pressure
options.drumFlag=0; %   [-] : if =1 then drum if =0 => no drum. 
options.eta_mec=0; %    [-] : mecanic efficiency of shafts bearings
         % options.comb is a structure containing combustion data : 
options.comb.Tmax=0; %      [�C] : maximum combustion temperature
options.comb.lambda=0; %    [-] : air excess
options.comb.x=0; %         [-] : the ratio O_x/C. Example 0.05 in CH_1.2O_0.05
options.comb.y=0; %         [-] : the ratio H_y/C. Example 1.2 in CH_1.2O_0.05
options.T_exhaust=0; %  [�C] : Temperature of exhaust gas out of the chimney
options.p_3=0; %        [-] : High pressure after last reheating
options.x4=0; %         [-] : Vapor ratio [gaseous/liquid] (in french : titre)
options.T_0=0; %        [�C] : Reference temperature
options.TpinchSub=0; %  [�C] : Temperature pinch at the subcooler
options.TpinchEx=0; %   [�C] : Temperature pinch at a heat exchanger
options.TpinchCond=0; % [�C] : Temperature pinch at condenser 
options.Tdrum=0; %      [�C] : minimal drum temperature
option.eta_SiC=0; %     [-] : Isotrenpic efficiency for compression
ption.eta_SiT=0; %     [-] : Isotrenpic efficiency for Turbine. It can be a vector of 2 values :
%             	             eta_SiT(1)=eta_SiT_HP,eta_SiT(2)=eta_SiT_others
display=0;
[ETA XMASSFLOW DATEN DATEX DAT MASSFLOW COMBUSTION FIG] = ST(P_e,options,display)