function [ETA, XMASSFLOW, DATEN, DATEX, DAT, MASSFLOW, COMBUSTION, FIG] = ST(P_e,options,display)
% ST Steam power plants modelisation
% ST(P_e,options,display) compute the thermodynamics states for a Steam
% power plant (combustion, exchanger, cycle) turbine based on several
% inputs (given in OPTION) and based on a given electricity production P_e.
% It returns the main results. It can as well plots graphs if input
% argument DISPLAY = true (<=> DISPLAY=1)
%
% INPUTS (some inputs can be dependent on others => only one of these 2 can
%         be activated)
% P_E = electrical power output target [kW]
% OPTIONS is a structure containing :
%   -options.nsout     [-] : Number of feed-heating
%   -options.reheat    [-] : Number of reheating
%   -options.T_max     [�C] : Maximum steam temperature
%   -options.T_cond_out[�C] : Condenseur cold outlet temperature
%   -options.p3_hp     [bar] : Maximum pressure
%   -options.drumFlag  [-] : if =1 then drum if =0 => no drum.
%   -options.eta_mec   [-] : mecanic efficiency of shafts bearings
%   -options.comb is a structure containing combustion data :
%       -comb.Tmax     [�C] : maximum combustion temperature
%       -comb.lambda   [-] : air excess
%       -comb.x        [-] : the ratio O_x/C. Example 0.05 in CH_1.2O_0.05
%       -comb.y        [-] : the ratio H_y/C. Example 1.2 in CH_1.2O_0.05
%   -options.T_exhaust [�C] : Temperature of exhaust gas out of the chimney
%   -options.p_3       [-] : High pressure after last reheating
%   -options.x4        [-] : Vapor ratio [gaseous/liquid] (in french : titre)
%   -options.T_0       [�C] : Reference temperature
%   -options.TpinchSub [�C] : Temperature pinch at the subcooler
%   -options.TpinchEx  [�C] : Temperature pinch at a heat exchanger
%   -options.TpinchCond[�C] : Temperature pinch at condenser
%   -options.Tdrum     [�C] : minimal drum temperature
%   -option.eta_SiC    [-] : Isotrenpic efficiency for compression
%   -option.eta_SiT    [-] : Isotrenpic efficiency for Turbine. It can be a vector of 2 values :
%             	             eta_SiT(1)=eta_SiT_HP,eta_SiT(2)=eta_SiT_others
% DISPLAY = 1 or 0. If 1, then the code should plot graphics. If 0, then
%          do not plot.
%
%OUPUTS :
% ETA is a vector with :
%   -eta(1) : eta_cyclen, cycle energy efficiency
%   -eta(2) : eta_toten, overall energy efficiency
%   -eta(3) : eta_cyclex, cycle exegy efficiency
%   -eta(4) : eta_totex, overall exergie efficiency
%   -eta(5) : eta_gen, Steam generator energy efficiency
%   -eta(6) : eta_gex, Steam generator exergy efficiency
%   -eta(7) : eta_combex, Combustion exergy efficiency
%   -eta(8) : eta_chemex, Chimney exergy efficiency (losses)
%   -eta(9) : eta_transex, Heat exchanger overall exergy efficiency
%   FYI : eta(i) \in [0;1] [-]
% Xmassflow is a vector with each feedheating massflow [kg/s] (respect to figure
%           2.33, page 91 "Thermal Power Plants" English version).
%           Xmassflow(1) = mass flow at 6_1 etc...
% DATEN is a vector with :
%   -daten(1) : perte_gen [kW]
%   -daten(2) : perte_mec [kW]
%   -daten(3) : perte_cond [kW]
% DATEX is a vector with :
%   -datex(1) : perte_mec    [kW]
%   -datex(2) : perte_totex  [kW]
%   -datex(3) : perte_rotex  [kW]
%   -datex(4) : perte_combex [kW]
%   -datex(5) : perte_condex [kW]
%   -datex(6) : perte_chemex [kW]
%   -datex(7) : perte_transex[kW]
% DAT is a matrix containing :
% dat = {T_1       , T_2       , ...       , T_6_I,     T_6_II, ... ;  [�C]
%        p_1       , p_2       , ...       , p_6_I,     p_6_II, ... ;  [bar]
%        h_1       , h_2       , ...       , h_6_I,     h_6_II, ... ;  [kJ/kg]
%        s_1       , s_2       , ...       , s_6_I,     s_6_II, ... ;  [kJ/kg/K]
%        e_1       , e_2       , ...       , e_6_I,     e_6_II, ... ;  [kJ/kg]
%        x_1       , x_2       , ...       , x_6_I,     x_6_II, ... ;   };[-]
% MASSFLOW is a vector containing :
%   -massflow(1) = m_a, air massflow [kg/s]
%   -massflow(2) = m_v, water massflow at 2 [kg/s]
%   -massflow(3) = m_c, combustible massflow [kg/s]
%   -massflow(4) = m_f, exhaust gas massflow [kg/s]
%
% COMBUSTION is a structure with :
%   -combustion.LHV    : the Lower Heat Value of the fuel [kJ/kg]
%   -combustion.e_c    : the combustible exergie         [kJ/kg]
%   -combustion.lambda : the air excess                   [-]
%   -combustion.Cp_g   : heat capacity of exhaust gas     [kJ/kg/K]
%   -combustion.fum  : is a vector of the exhaust gas composition :
%       -fum(1) = m_O2f  : massflow of O2 in exhaust gas [kg/s]
%       -fum(2) = m_N2f  : massflow of N2 in exhaust gas [kg/s]
%       -fum(3) = m_CO2f : massflow of CO2 in exhaust gas [kg/s]
%       -fum(4) = m_H2Of : massflow of H2O in exhaust gas [kg/s]
%
% FIG is a vector of all the figure you plot. Before each figure, define a
% figure environment such as:
%  "FIG(1) = figure;
%  plot(x,y1);
%  [...]
%   FIG(2) = figure;
%  plot(x,y2);
%  [...]"
%  Your vector FIG will contain all the figure plot during the run of this
%  code (whatever the size of FIG).

%% YOUR WORK

% Exemple of how to use 'nargin' to check your number of inputs
if nargin<3
    display = 1;
    if nargin<2
        options = struct();
        if nargin<1
            P_e = 288e3; % [kW] Puissance �nerg�tique de l'installation
        end
        options.nsout=8; %   [-] : Number of feed-heating %%
        options.reheat=1; %    [-] : Number of reheating %%
        options.T_max=565; %     [�C] : Maximum steam temperature %%
        options.T_cond_out=24; %[�C] : Condenseur cold outlet temperature %%
        options.p3_hp=310; %     [bar] : Maximum pressure  %%
        options.drumFlag=1; %   [-] : if =1 then drum if =0 => no drum.  %%
        options.eta_mec=0.99; %    [-] : mecanic efficiency of shafts bearings %%
        % options.comb is a structure containing combustion data :
        options.comb=struct;
        options.comb.Tmax=1050; %      [�C] : maximum combustion temperature
        %options.comb.lambda=0; %    [-] : air excess
        options.comb.x=0; %         [-] : the ratio O_x/C. Example 0.05 in CH_1.2O_0.05
        options.comb.y=0; %         [-] : the ratio H_y/C. Example 1.2 in CH_1.2O_0.05
        options.T_exhaust=80; %  [�C] : Temperature of exhaust gas out of the chimney %%
        options.p_3=62; %        [-] : High pressure after last reheating  %%
        options.x4=0.99; %         [-] : Vapor ratio [gaseous/liquid] (in french : titre) %%
        options.T_0=15; %        [�C] : Reference temperature %%
        options.TpinchSub=4; %  [�C] : Temperature pinch at the subcooler %%
        options.TpinchEx=10; %   [�C] : Temperature pinch at a heat exchanger %%
        options.TpinchCond=8; % [�C] : Temperature pinch at condenser %%
        options.Tdrum=148.7; %      [�C] : minimal drum temperature  %%
        options.eta_SiC=0.85; %     [-] : Isotrenpic efficiency for compression %%
        options.eta_SiT=[0.89 0.89]; %     [-] : Isotrenpic efficiency for Turbine. It can be a vector of 2 values :
        %             	             eta_SiT(1)=eta_SiT_HP,eta_SiT(2)=eta_SiT_others
        %             	             %%
    end
end


% Exemple of how to use (isfield' to check if an option has been given (or
% not)
if isfield(options,'T_0')
    T_0 = options.T_0+273.15;
else
    T_0 = 288.15;  % [K]
end

if ~isfield(options,'drumFlag')
    options.drumFlag=1;
end

if ~isfield(options,'eta_mec')
    options.eta_mec=0.99;
end

if ~isfield(options,'nsout')
    options.nsout=8;
end

if ~isfield(options,'reheat')
    options.reheat=1;
end

if ~isfield(options,'T_max')
    options.T_max=565; % [�C]
end

if ~isfield(options,'T_cond_out')
    options.T_cond_out=26; % [�C]
end

if ~isfield(options,'p3_hp')
    options.p3_hp=310; % [bar]
end

if ~isfield(options,'comb')
    options.comb=struct;
    options.comb.Tmax=1400; %      [�C] : maximum combustion temperature
    %options.comb.lambda=0; %    [-] : air excess
    options.comb.x=0; %         [-] : the ratio O_x/C. Example 0.05 in CH_1.2O_0.05
    options.comb.y=4; %         [-] : the ratio H_y/C. Example 1.2 in CH_1.2O_0.05
else
    if ~isfield(options.comb,'Tmax')
        options.Tmax=565; % [�C]
    end
    if ~isfield(options.comb,'x')
        options.x=0;
    end
    if ~isfield(options.comb,'y')
        options.y=0;
    end
end

if ~isfield(options,'T_exhaust')
    options.T_exhaust=80; % [�C]
end

if ~isfield(options,'eta_SiC')
    options.eta_SiC=0.85;
end

if ~isfield(options,'eta_SiT')
    options.eta_SiT=[0.89 0.89];
end

if (~isfield(options,'p_3') && isfield(options,'x4'))
    options.p_3=detente(eta_SiT(1),options.TpinchCond+options.T_cond_out,...
        options.x4,options.T_max); % [bar]
else
    options.p_3=0.2*options.p3_hp;
end

if ~isfield(options,'x4')
    options.x4=0.99; % [-]
end

if ~isfield(options,'p_3')
    options.p_3=7e6; % [bar]
end

if ~isfield(options,'TpinchSub')
    options.TpinchSub=4; % [�C]
end

if ~isfield(options,'TpinchEx')
    options.TpinchEx=8; % [�C]
end

if ~isfield(options,'TpinchCond')
    options.TpinchCond=10; % [�C]
end

if ~isfield(options,'Tdrum')
    options.Tdrum=148.7; % [�C]
end

if ~isfield(options,'kcc')
    options.k_cc=6.2/7; % [-] : Coefficient of pressure losses in the boiler
end

if ~isfield(options,'GraphCompute')
    options.GraphCompute=1;
end

if ~isfield(options,'erCompute')
    options.erCompute=0;
end

%
nsout=options.nsout;
eta_SiT=options.eta_SiT;
TpinchSub=options.TpinchSub;
TpinchEx=options.TpinchEx;
eta_SiC=options.eta_SiC;
k_cc= options.k_cc; % [-] : Coefficient of pressure losses in the boiler

%
%Prealocations
t6=zeros(nsout+1,1);
p6=zeros(nsout+1,1);
h6=zeros(nsout+1,1);
s6=zeros(nsout+1,1);
x6=zeros(nsout+1,1);
x6=zeros(nsout+1,1);

h0=XSteam('hL_T',T_0-273.15);
s0=XSteam('sL_T',T_0-273.15);

t3=options.T_max; % [�C]
p3=options.p3_hp; % [bar]
h3=XSteam('h_pT',p3,t3);
s3=XSteam('s_pT',p3,t3);
x3=1; % Vapeur surchauf�e
e3=(h3-h0)-T_0*(s3-s0); % [kJ/kg]

%Etat 4
if options.reheat==1
    p4=options.p_3/k_cc; % Presure drop in the boiler
else
    p4=options.p_3;
end
[ti,pi,hi,si,xi] = expansion(eta_SiT(1),t3,p3,h3,s3,x3,p4,0);
t4=ti(end);
p4=pi(end);
h4=hi(end);
s4=si(end);
x4=xi(end);
e4=(h4-h0)-T_0*(s4-s0);

%Etat 5
if options.reheat==1
    t5=t3;
    p5=options.p_3;
    h5=XSteam('h_pT',p5,t5);
    s5=XSteam('s_pT',p5,t5);
    e5=(h5-h0)-T_0*(s5-s0);
    x5=1; % vapeur surchaf�e
else % Etat 4 = Etat 5
    t5=t4;
    p5=p3;
    h5=h4;
    s5=s4;
    x5=x4;
    e5=e4;
end

%Etat 7
t7=options.T_cond_out+options.TpinchCond;
p7=XSteam('psat_T',t7);
h7=XSteam('hL_T',t7);
s7=XSteam('sL_T',t7);
e7=(h7-h0)-T_0*(s7-s0);
x7=0;

%Etat(s) 6
[t_i,p_i,h_i,s_i,x_i] = expansion(eta_SiT(end),t5,p5,h5,s5,x5,p7,nsout-1);
t6=[ flip(t_i(2:end)) t4];
p6=[ flip(p_i(2:end)) p4];
h6=[ flip(h_i(2:end)) h4];
s6=[ flip(s_i(2:end)) s4];
x6=[ flip(x_i(2:end)) x4];
e6=(h6-h0)-T_0*(s6-s0);

if options.drumFlag==1
    ind_drum=1;
    t_sat=0;
    while t_sat < options.Tdrum % Trouve le soutirage allant jusqu'au drum
        ind_drum=ind_drum+1;
        t_sat=XSteam('Tsat_p',p6(ind_drum));
    end
    p8=XSteam('psat_T',options.Tdrum);
else  % If there's no drum
    NPSH_no_drum=4; % [m]
    secu_marg=1.5; % Margin factor of the NPSH
    rho_eau=1e3; % [kg/m^3]
    p8=p6(end)+ rho_eau*9.81*(NPSH_no_drum+secu_marg)*1e-5;
    ind_drum=0;
end

%Etat 8
[t8,h8,s8,x8,e8] = compression(options.eta_SiC,p8,h7(1),s7(1));

% Etat(s) 7 et 9

if options.drumFlag==1 % If there's a drum, we focus on the pre-drum feed_heatings
    p_feed=p6(2:ind_drum-1);
    h_feed=h6(2:ind_drum-1);
else
    p_feed=p6(2:end);
    h_feed=h6(2:end);
end

t90=40; % [�C]  we pose a post-subcooler temperature

%Before the drum
[t_fs,h_fs,p_fs,s_fs,t9,h9,p9,s9] = states7_9_pre_drum2(h_feed,p_feed,t90);

if options.drumFlag==1
    p_drum=XSteam('psat_T',options.Tdrum);
    t7=[t7 ; t_fs ; options.Tdrum];
    p7=[p7 ; p_fs ; p_drum];
    h7=[h7; h_fs ; XSteam('hL_p',p_drum)];
    s7=[s7 ; s_fs ; XSteam('sL_p',p_drum)];
else
    t7=[t7 ; t_fs ];
    p7=[p7 ; p_fs ];
    h7=[h7; h_fs ];
    s7=[s7 ; s_fs];
end
% After the drum if there's one
if options.drumFlag==1
    p_feed=p6(ind_drum+1:end);
    h_feed=h6(ind_drum+1:end);
    [t_w,h_w,p_w,s_w,t_9,h_9,p_9,s_9] = states7_9_post_drum2(h_feed,p_feed,h7(end),s7(end),options.eta_SiC);
    
    t7=[t7 ; t_w ];
    p7=[p7 ; p_w ];
    h7=[h7; h_w];
    s7=[s7 ; s_w];
    
    t9=[t9 ; t_9 ];
    p9=[p9 ; p_9 ];
    h9=[h9; h_9];
    s9=[s9 ; s_9];
end
e7= (h7-h0)-T_0*(s7-s0);
e9= (h9-h0)-T_0*(s9-s0);
x7=zeros(size(t7));
x9=zeros(size(t9));

%Etat 1
t1=t9(end);
p1=p9(end);
s1=s9(end);
h1=h9(end);
e1=e9(end);
x1=0;

% Etat 2
p2=p3/k_cc;
[t2,h2,s2,x2,e2] = compression(eta_SiC,p2,h1,s1);

% Flowrates X
if options.nsout >0
[X6_flows,h_sc,e_sc]=X_flowrates();
else
    X6_flows=0;
end

m6 = flow_mass(X6_flows);

m7=(1+sum(X6_flows))*m6;
m5=(1+sum(X6_flows(1:end-1)))*m6;
m2=m7;
m3=m2;

XMASSFLOW=m6*X6_flows;

%Prealocations
DATEN=zeros(3,1);
DATEX=zeros(7,1);
ETA=zeros(9,1);
MASSFLOW=zeros(4,1);

% Combustion
[COMBUSTION,h_fum,ratio_air,ratio_fum,ma1,h_air,TmaxComb] = combustion(options.comb);
Q_I=m3*(h3-h2)+ m5*(h5-h4); % Combustion heat energy
h_exh = getEnthalpy('fum',ratio_fum,T_0,options.T_exhaust+273.15);
mc=Q_I/((1+COMBUSTION.lambda*ma1)*(h_fum-h_exh));
m_air = COMBUSTION.lambda*ma1*mc;
mf=m_air+mc;

m2=m3; % flow rate a point 2

COMBUSTION.fum=ratio_fum.*mf;
ef = h_fum - T_0*(integral(@(T) getCpGas('fum',ratio_fum,T)./T,T_0,TmaxComb));
e_exh = h_exh - T_0*(integral(@(T) getCpGas('fum',ratio_fum,T)./T,T_0,options.T_exhaust+273.15));

if options.erCompute==1 % Preheating of the air in the boiler
    lambd=COMBUSTION.lambda;
    Tb=t2+TpinchEx+273.15;
    Tc=options.T_exhaust;
    Tin=273.15+15; % Air enters at 15�C
    
    eq=@(T) (m_air*getEnthalpy('air',ratio_air,Tin,T))-...
        (mf*getEnthalpy('fum',ratio_fum,Tc,Tb));
    
    optis = optimset('Display','off');
    Ta1=fsolve(eq,300,optis);
    
    ea=getEnthalpy('air',ratio_air,Tin,Ta1)-T_0*integral(@(T) getCpGas('air',ratio_air,T)./T,T_0,Ta1);
    er=ea*(lambd*ma1)/(lambd*ma1+1);
else % no preheating of the air in the boiler
    er=0;
end

% Puissances
P_hp=m3*(h3-h4); % HP turbine
P_hpx=m3*(e3-e4);
P_ip=m6*(1+sum(X6_flows(1:ind_drum-2)))*(h5-h6(ind_drum-1)); % IP turbine
P_ipx=m6*(1+sum(X6_flows(1:ind_drum-2)))*(e5-e6(ind_drum-1));
for k=ind_drum:nsout-1
   P_ip= P_ip  +X6_flows(k-1)*m6*(h5-h6(k)); 
   P_ipx=P_ipx +X6_flows(k-1)*m6*(e5-e6(k));  
end
P_lp=m6*(h6(ind_drum)-h6(1));% LP Turbine
P_lpx=m6*(e6(ind_drum)-e6(1));
for k=2:ind_drum-1
   P_lp= P_lp  +X6_flows(k-1)*m6*(h6(ind_drum)-h6(k)); 
   P_lpx=P_lpx +X6_flows(k-1)*m6*(e6(ind_drum)-e6(k)); 
end

Ppa=m2*(h2-h1); % Pump Pa
Ppax=m2*(e2-e1); 
Ppe=m6*(h8-h7(1)); % Pump Pe
Ppex=m6*(e8-e7(1));

if options.drumFlag==1
    Ppb= m2*(h9(ind_drum)-h7(ind_drum));
    Ppbx= m2*(e9(ind_drum)-e7(ind_drum));    
    DATEN(3)=DATEN(3)+sum(XMASSFLOW(1:ind_drum-2)*(h_sc-h7(1))); % perte_cond [kW]
else
    DATEN(3)=DATEN(3); % perte_cond [kW]
    Ppb=0;
    Ppbx=0;
end

PmC=Ppe+Ppa+Ppb;
PmCx=Ppex+Ppax+Ppbx;
%PmT=P_e/options.eta_mec;
PmT=P_hp+P_ip+P_lp;
PmTx=P_hpx+P_ipx+P_lpx;

Q_Ix=m3*(e3-e2)+ m5*(e5-e4); % exergy of the combustion

%Massflows
MASSFLOW(1)=m_air;
MASSFLOW(2)=m2;
MASSFLOW(3)=mc;
MASSFLOW(4)=mf;

% Pertes
DATEN(1)=mc*LHV-(m3*(h3-h2)+m5*(h5-h4));  % perte_gen [kW]
DATEN(2)=P_e*(1/options.eta_mec - 1); % perte_mec [kW]
DATEN(3)=m6*(h6(1)-h7(1)); % perte_cond [kW]

DATEX(1)=P_e*(1/options.eta_mec - 1); % perte_mec [kW]
DATEX(2)=mc*ec-P_e; % perte_totex [kW]
DATEX(3)= (PmTx-PmCx) - (PmT-PmC); % perte_rotex [kW]
DATEX(4)=mc*ec - mf*(ef-er); % perte_combex [kW]
DATEX(5)= m6*(e6(1)-e7(1)); % perte_condex [kW]
DATEX(6)=mf*(e_exh - er); % perte_chemex [kW]
DATEX(7)=mf*(ef-e_exh) - m3*(e3-e2);  % perte_transex [kW]

if options.drumFlag==1
    DATEX(5)=DATEX(5)+sum(XMASSFLOW(1:ind_drum-1)*e_sc-e7(1)); % perte_condex [kW]
else
    DATEX(5)=DATEX(5); % perte_condex [kW]
end

MASSFLOW(2)=m2; % water massflow at 2 [kg/s]
MASSFLOW(3)=mc; % combustible massflow [kg/s]
MASSFLOW(4)=mf; % exhaust gas massflow [kg/s]

DAT=[ t1 t2 t3 t4 t5 t6 t7' t8 t9' ;...
    p1 p2 p3 p4 p5 p6 p7' p8 p9' ;...
    h1 h2 h3 h4 h5 h6 h7' h8 h9' ;...
    s1 s2 s3 s4 s5 s6 s7' s8 s9' ;...
    e1 e2 e3 e4 e5 e6 e7' e8 e9' ;...
    x1 x2 x3 x4 x5 x6 x7' x8 x9'];

ETA(1)= (PmT-PmC)/Q_I; % eta_cyclen
ETA(2)= P_e/(mc*LHV); % eta_toten
ETA(3)= (PmT-PmC)/Q_Ix; % eta_cyclex
ETA(4)=  P_e/(mc*ec); % eta_totex
ETA(5)= Q_I/(mc*LHV); % eta_gen
ETA(6)= Q_Ix/(mc*ec); % eta_gex
ETA(7)= mf*(ef-er)/(mc*ec); % eta_combex
ETA(8)= (ef - e_exh)/(ef-er); % eta_chemex
ETA(9)= Q_Ix/(mf*(ef-e_exh)); %eta_transex

% Plots
if options.GraphCompute==1
if display ==1
    disp='on';
else
    disp='off';
end
%Cloche T-S
S_plot=linspace(0,9.156,300);
T_HS=zeros(size(S_plot));
H_HS=zeros(size(S_plot));
for c=1:length(T_HS)
    T_HS(c)=XSteam('Tsat_s',S_plot(c));
end
[~,ind_tmax]=max(T_HS);
for c=1:ind_tmax
    H_HS(c)=XSteam('hL_T',T_HS(c));
end
for c=ind_tmax+1:length(T_HS)
    H_HS(c)=XSteam('hV_T',T_HS(c));
end

% cells (It's easier to plot with cells)
s_cell= { s1; s2 ; s3 ; s4 ; s5 ; s6' ; s7 ;  s8 ; s9 };
t_cell= { t1; t2 ; t3 ; t4 ; t5 ; t6' ; t7 ;  t8 ; t9 };
h_cell= { h1; h2 ; h3 ; h4 ; h5 ; h6' ; h7 ;  h8 ; h9 };

FIG = gobjects(4,1);

% T-S Diagram
FIG(1)=figure('visible',disp);
hold on
plot(S_plot,T_HS);
leg=cell(10,1);
leg{1}='Saturation curve';
for c=1:9
    plot(s_cell{c},t_cell{c},'*')
    leg{c+1}=['State ' num2str(c)];
end
legend(leg)
title('Steam turbine T-s diagram');
ylabel('Temperature [�C]');
xlabel('Entropy [kJ/kg/K]');

for c=1:9
    plot_TS(c);
end

% H-S Diagram
FIG(2)=figure('visible',disp);
hold on
plot(S_plot,H_HS);
title('Steam turbine h-s diagram');
ylabel('Enthalpie [kJ/kg]');
xlabel('Entropy [kJ/kg/K]');
legh=cell(10,1);
legh{1}='Saturation curve';
for c=1:9
    plot(s_cell{c},h_cell{c},'*')
    legh{c+1}=['State ' num2str(c)];
end
legend(legh)

for c=1:9
    plot_HS(c);
end

FIG(3)=figure('visible',disp);
pie([P_e ; DATEN],{sprintf('%s \n %.1f [MW]','Effective power',P_e*1e-3),...
    sprintf('%s \n %.1f [MW]','Generator losses',DATEN(1)*1e-3),...
    sprintf('%s \n %.1f [MW]','Mechanical losses',DATEN(2)*1e-3),...
    sprintf('%s \n %.1f [MW]','Condensor loss',DATEN(3)*1e-3)});
title(sprintf('Steam Turbine Energy pie chart \n Primary power :  %.1f [MW]',mc*LHV*1e-3));


FIG(4)=figure('visible',disp);
pie([P_e ; DATEX(1); DATEX(3:end)],{sprintf('%s \n %.1f [MW]','Effective power',P_e*1e-3),...
    sprintf('%s \n %.1f [MW]','Mechanical losses',DATEX(1)*1e-3),...
    sprintf('%s \n %.1f [MW]','Turbine & compressor irreversibilities',DATEX(3)*1e-3),...
    sprintf('%s \n %.1f [MW]','Combustion irreversibility',DATEX(4)*1e-3)...
    sprintf('%s \n %.1f [MW]','Condensor loss',DATEX(5)*1e-3),...
    sprintf('%s \n %.1f [MW]','Chemney loss',DATEX(6)*1e-3),sprintf('%s : %.1f [MW]','Heat transfer loss',DATEX(7)*1e-3)});
title(sprintf('Steam Turbine Exergy pie chart \n Primary exergy flux :  %.1f [MW]',mc*ec*1e-3));

end

    function [ti,pi,hi,si,xi] = expansion(eta_SiT,t_in,p_in,h_in,s_in,x_in,p_out,nsout)
        % ti,pi,hi,si,xi: state of the feed-heatings
        % t_in,p_in,h_in,s_in : conditions at the turbine inlet
        % p_out pressure after expansion
        % nsout : number of feedheatings
        
        % "main" expansion
        s_outs=s_in;
        h_outs=XSteam('h_ps',p_out,s_outs); % Isentropic expansion
        h_out= h_in + eta_SiT*(h_outs - h_in); % Real expansion
        t_out= XSteam('t_ph',p_out,h_out);
        s_out= XSteam('s_ph',p_out,h_out);
        x_out= XSteam('x_ph',p_out,h_out);
        
        % Feed-heatings
        step_h=(h_out-h_in)/(nsout+1);
        hi=h_in:step_h:h_out; % Enthalpies of the differents feedheatings
        % Prealocations
        ti=zeros(size(hi));
        pi=zeros(size(hi));
        si=zeros(size(hi));
        xi=zeros(size(hi));
        
        ti(1)=t_in;
        ti(end)=t_out;
        pi(1)=p_in;
        pi(end)=p_out;
        xi(1)=x_in;
        xi(end)=x_out;
        si(1)=s_in;
        si(end)=s_out;
        
        if nsout > 0
            for i=1:nsout
                [ti(i+1),pi(i+1),si(i+1),xi(i+1)]=expan_inter(eta_SiT,hi(i+1),hi(1),si(1));
            end
        end
    end

    function [t2,h2,s2,x2,e2] = compression(eta_SiC,p2,h1,s1)
        % Adiabatic compression
        
        h_2s=XSteam('h_ps',p2,s1);
        h2=h1+1/eta_SiC*(h_2s-h1);
        s2=XSteam('s_ph',p2,h2);
        x2=XSteam('x_ph',p2,h2);
        t2=XSteam('T_ph',p2,h2);
        e2=(h2-h0)-T_0*(s2-s0);
    end

    function [t_out,p_out,s_out,x_out]=expan_inter(eta_SiT,h_out,h_init,s_i)
        % Compute the intermediate states of the feed-heatings
        
        h_outs=(h_out-(1-eta_SiT)*h_init)/(eta_SiT); % ISentropic expansion
        s_outs=s_i;
        p_out=XSteam('p_hs',h_outs,s_outs);
        t_out=XSteam('t_ph',p_out,h_out);
        s_out=XSteam('s_ph',p_out,h_out);
        x_out = XSteam('x_ps',p_out,s_out);
        
        if x_out < 0.88
            error('The turbine can only work with x > 0.88')
        end
    end

    function [t_w,h_w,p_w,s_w,t_9,h_9,p_9,s_9] = states7_9_pre_drum2(h6,p6,t90)
        %Computes the datas of states 7 and 9 and the rates X
        
        ns=length(h6); % Number of feed-heating used to warm up the water
        %Prealocations
        t_9=zeros(ns+1,1);
        p_9=zeros(ns+1,1);
        s_9=zeros(ns+1,1);
        h_9=zeros(ns+1,1);
        t_w=zeros(ns,1);
        p_w=zeros(ns,1);
        s_w=zeros(ns,1);
        h_w=zeros(ns,1);
        
        p_9(1:end)=p8; % presure is the same
        
        for i=1:ns
            % State 7
            p_w(i)=p6(i);  % The heat exchange is isobaric
            t_w(i)=XSteam('tsat_p',p_w(i)); % After the heat exchange, the water is saturated
            h_w(i)=XSteam('hL_p',p_w(i));
            s_w(i)=XSteam('sL_p',p_w(i));
            %State 9
            t_9(i+1)=t_w(i)-TpinchEx;
            h_9(i+1)=XSteam('h_pT',p_9(i+1),t_9(i+1));
            s_9(i+1)=XSteam('s_pT',p_9(i+1),t_9(i+1));
        end
        t_9(1)=t90;
        h_9(1)=XSteam('h_pT',p_9(1),t_9(1));
        s_9(1)=XSteam('s_pT',p_9(1),t_9(1));
        
    end

    function [t_w,h_w,p_w,s_w,t_9,h_9,p_9,s_9] = states7_9_post_drum2(h6,p6,h_bef,s_bef,eta_SiC)
        % Computes the datas of states 7 and 9 and the rates X after the
        % drum
        
        ns=length(h6); % Number of feed-heating used to warm up the water
        %Prealocations
        t_9=zeros(ns+1,1);
        p_9=zeros(ns+1,1);
        s_9=zeros(ns+1,1);
        h_9=zeros(ns+1,1);
        t_w=zeros(ns,1);
        p_w=zeros(ns,1);
        s_w=zeros(ns,1);
        h_w=zeros(ns,1);
        
        for i=1:ns
            % State 7
            p_w(i)=p6(i);
            t_w(i)=XSteam('tsat_p',p_w(i));
            h_w(i)=XSteam('hL_p',p_w(i));
            s_w(i)=XSteam('sL_p',p_w(i));
            
            t_9(i+1)=t_w(i)- max(TpinchEx,TpinchSub); %State 9
        end
        
        NPSH=4; %[m]
        secu_margin=1.5; % Margin factor of the NPSH
        t_Pa=t_9(end); % [�C] Temperature at the pump Pa
        rho=1e3; % [kg/m^3]
        p_sat=XSteam('psat_T',t_Pa);
        g=9.81; % [m/s^2]
        p_9(1:end)=p_sat+ rho*g*(NPSH+secu_margin)*1e-5;
        
        for i=2:ns+1
            h_9(i)=XSteam('h_pT',p_9(i),t_9(i));
            s_9(i)=XSteam('s_pT',p_9(i),t_9(i));
        end
        
        % Compression at Pb
        [t_9(1),h_9(1),s_9(1),~,~] = compression(eta_SiC,p_9(1),h_bef,s_bef);
        
        
        % Computations of the X
        % They can be computed by resolving a linear system of type Ax=B
        % Here x= [ X_i X_i+1 ... X_N]
        
    end

    function [X,hsc,esc]= X_flowrates()      
        
        % State sc : Between the subcooler and the Condensor
        tsc=t8+TpinchSub;
        psc=p9(1);
        hsc=XSteam('h_pt',psc,tsc);
        ssc=XSteam('s_pt',psc,tsc);
        esc=(hsc - h0) - T_0*(ssc-s0);
        
        % Computations of the X
        % They can be compbuted by resolving a linear system of type Ax=B
        % Here x= [X_I X_II ...]
        
        A=zeros(nsout,nsout);
        B=zeros(nsout,1);
        
        if options.drumFlag==1
           id=ind_drum-1; % drum indice
           
           %Before the drum
           for k=1:id-1
           kd=k+1;
           % Before the drum
           A(k,k)= h6(kd)-h7(kd);
           A(k,k+1:id-1)=h7(kd+1)-h7(kd);
           A(k,1:id-1)=A(k,1:id-1) - (h9(kd)-h9(kd-1));
           
           B(k)=h9(kd)-h9(kd-1);
           end
           
           %At the drum
           k=id;
           kd=k+1;
           A(k,k)= h7(kd)-h6(kd);
           A(k,k+1:nsout)= A(k,k+1:nsout)+(h7(kd)-h7(kd+1));
           A(k,1:k-1)= A(k,1:k-1) - (h7(kd)-h9(kd-1)); 
           
           B(k)=h9(kd-1)- h7(kd);
           
           % After the drum
           
           for k=id+1:nsout-1
               kd=k+1;
               % Before the drum
               A(k,k)= h6(kd)-h7(kd);
               A(k,k+1:nsout)=  A(k,k+1:nsout) + (h7(kd+1)-h7(kd));
               A(k,1:nsout)=A(k,1:nsout) - (h9(kd)-h9(kd-1));
               
               B(k)=h9(kd)-h9(kd-1);
           end
           
           % Last exchanger
           k=nsout;
           kd=k+1;
           A(k,k)= h6(kd)-h7(kd);
           A(k,1:nsout)= A(k,1:nsout) - (h9(kd)-h9(kd-1));
             
            B(k)=h9(kd)-h9(kd-1);
            
            X=A\B;
            
        else
            
        end
    end

    function m6 = flow_mass(X)
        % Computes the flow rate at point 7
        
        % Normalisation of Wm with regard to the point 7
        %Hp Turbine and start of the IP
        if options.nsout > 0
        Wm=(1+sum(X))*(h3-h4)+ (1+sum(X(1:end-1)))*(h5-h6(end-1));
        tmp_deb=(1+sum(X(1:end-1)));
        else
           Wm=(h3-h4)+(h5-h6);
        end
            for k=nsout-1:-1:2
            tmp_deb=tmp_deb -X(k);
            Wm=Wm+tmp_deb*(h6(k)-h6(k-1));
            end       
         m6=P_e/(options.eta_mec*Wm); % Flow rate a point 7
    end

    function [COMBUSTION,h_fum,ratio_air,ratio_fum,ma1,h_air,Tmax] = combustion(comb)
        x=comb.x;
        y=comb.y;
        ma1=(1+(y-2*x)/4)*(32+3.76*28.15)/(12.01+1.008*y+16*x);
        LHV=getLHV(x,y); % [kJ/kg]
        cp_fuel=getCpfuel(x,y); % [kJ/kg/K]
        ec=getFuelExerfy(x,y);% [kJ/kg]
        T_ext=273.15+ 15; % [K]
        h_fuel=cp_fuel * (T_ext-T_0);
        
        Mm_air=0.21*32+0.79*28; % [g/mol]
        ratio_O2=0.21*32/Mm_air;
        ratio_N2=0.79*28/Mm_air;
        ratio_air=[ratio_O2 ; ratio_N2];
        h_air = getEnthalpy('air',ratio_air,T_0,T_ext);
        
        opts=optimoptions('fsolve','Display','off');
        if isfield(comb,'lambda')
            lambda = comb.lambda;
            Tmax = fsolve(@(Tmax) eqLambdaTmax(lambda,x,y,ma1,Tmax,LHV,h_air,h_fuel),1300,opts);
        elseif isfield(comb,'Tmax')
            Tmax = comb.Tmax + 273.15;
            lambda = fsolve(@(lambda) eqLambdaTmax(lambda,x,y,ma1,Tmax,LHV,h_air,h_fuel),1.5,opts);
        else
            error('comb must contain eitheir lambda or Tmax')
        end
        
        Mm_fum=(12+y+16*x+lambda*(1+0.25*(y-2*x))*(32+3.76*28));
        ratio_O2f=((lambda-1)*(1+0.25*(y-2*x)))*32/Mm_fum;
        ratio_N2f=(1+0.25*(y-2*x))*3.76*lambda*28/Mm_fum;
        ratio_CO2f=44/Mm_fum;
        ratio_H2Of=0.5*y*18/Mm_fum;
        ratio_fum=[ratio_O2f; ratio_N2f; ratio_CO2f; ratio_H2Of];
        h_fum = getEnthalpy('fum',ratio_fum,T_0,Tmax);
        cpg = getCpGas('fum',ratio_fum,400);
        
        COMBUSTION=struct;
        COMBUSTION.LHV=LHV;
        COMBUSTION.e_c=ec;
        COMBUSTION.lambda=lambda;
        COMBUSTION.cp_g=cpg;
    end

    function LHV =getLHV(x,y)
        LHV=(393400+102250*y-(x/(1+0.5*y))*(111000+102250*y))/(12+y+x*16);
    end

    function cpFuel=getCpfuel(x,y)
        Mm=12+y+16*x; % [kg/kmole]
        if x==0
            switch y
                case 0 % C
                    cpFuel=10.4/Mm;
                case 1.8 % CH1.8
                    cpFuel=14.4/Mm;
                case 4 % CH4
                    cpFuel=35.3/Mm;
            end
        elseif (y==0)&&(x==1) % CO
            cpFuel=29.1/Mm;
        else % CH4
            cpFuel=35.3/Mm;
        end
    end

    function ec=getFuelExerfy(x,y)
        if x==0
            switch y
                case 0 % C
                    ec=34160; % [kJ/kg]
                case 1.8 % CH1.8
                    ec=45710; % [kJ/kg]
                case 4 % CH4
                    ec=52215; % [kJ/kg]
            end
        elseif (y==0)&&(x==1) % CO
            ec=9845; % [kJ/kg]
        else % CH4
            ec=52215; % [kJ/kg]
        end
    end

    function cp = getCpGas(gas,composition,T)
        switch gas
            case 'air'
                if T>=300
                    cp=composition(1)*janaf('c','O2',T)+composition(2)*janaf('c','N2',T);
                elseif T<=300
                    cp=ones(size(T))*(composition(1)*janaf('c','O2',300)+...
                        composition(2)*janaf('c','N2',300));
                else
                    ind=find(T<300,1,'last');
                    cp=[(ones(1,ind)*(composition(1)*janaf('c','O2',300)+...
                        composition(2)*janaf('c','N2',300)))...
                        (composition(1)*janaf('c','O2',T(ind+1:end))+...
                        composition(2)*janaf('c','N2',T(ind+1:end)))];
                end
            case 'fum'
                if T>=300
                    cp=composition(1)*janaf('c','O2',T)+...
                        composition(2)*janaf('c','N2',T)+...
                        composition(3)*janaf('c','CO2',T)+...
                        composition(4)*janaf('c','H2O',T);
                elseif T <= 300
                    cp=ones(size(T))*(composition(1)*janaf('c','O2',300)+...
                        composition(2)*janaf('c','N2',300)+...
                        composition(3)*janaf('c','CO2',300)+...
                        composition(4)*janaf('c','H2O',300));
                else
                    ind = find(T < 300,1,'last');
                    cp = [(ones(1,ind)*(composition(1)*janaf('c','O2',300)+...
                        composition(2)*janaf('c','N2',300)+...
                        composition(3)*janaf('c','CO2',300)+...
                        composition(4)*janaf('c','H2O',300))),...
                        composition(1)*janaf('c','O2',T(ind+1:end))+...
                        composition(2)*janaf('c','N2',T(ind+1:end))+...
                        composition(3)*janaf('c','CO2',T(ind+1:end))+...
                        composition(4)*janaf('c','H2O',T(ind+1:end))];
                    
                end
            otherwise
                error('Gas must be air ou flue gass');
        end
    end

    function h = getEnthalpy(type,compo,T1,T2)
        h = integral(@(T) getCpGas(type,compo,T),T1,T2);
    end

    function eq =eqLambdaTmax(lambda,x,y,ma1,Tmax,LHV,h_air,h_fuel)
        Mm_fum=(12+y+16*x+lambda*(1+0.25*(y-2*x))*(32+3.76*28));
        ratio_O2f=((lambda-1)*(1+0.25*(y-2*x)))*32/Mm_fum;
        ratio_N2f=(1+0.25*(y-2*x))*3.76*lambda*28/Mm_fum;
        ratio_CO2f=44/Mm_fum;
        ratio_H2Of=0.5*y*18/Mm_fum;
        
        compo_fum = [ratio_O2f; ratio_N2f; ratio_CO2f; ratio_H2Of];
        hfum = getEnthalpy('fum',compo_fum,T_0,Tmax);
        eq=(lambda*ma1+1)*hfum-lambda*ma1*h_air-LHV-h_fuel;
    end

    function plot_TS(ind)
        n_pt=100;
        switch ind
            case 1 % Line 1-2
                p=linspace(p1,p2,n_pt);
                t=zeros(size(p));
                s=zeros(size(p));
                s(1)=s1;
                t(1)=t1;
                for k=2:n_pt-1
                    [t(k),~,s(k),~,~] = compression(eta_SiC,p(k),h1,s1);
                end
                t(end)=t2;
                s(end)=s2;
                plot(s,t,'k','HandleVisibility','off');
            case 2 % Line 2-3
                p=linspace(p2,p3,n_pt);
                s=linspace(s2,s3,n_pt);
                t=zeros(size(s));
                t(1)=t2;
                for k=2:n_pt-1
                    t(k)=XSteam('T_ps',p(k),s(k));
                end
                t(end)=t3;
                plot(s,t,'k','HandleVisibility','off');
            case 3 % Line 3-4
                p=linspace(p3,p4,n_pt);
                t=zeros(size(p));
                s=zeros(size(p));
                t(1)=t3;
                s(1)=s3;
                for k=2:n_pt-1
                    [T4,~,~,S4,~] = expansion(eta_SiT(1),...
                        t3,p3,h3,s3,x3,p(k),0);
                    t(k)=T4(end);
                    s(k)=S4(end);
                end
                t(end)=t4;
                s(end)=s4;
                plot(s,t,'k','HandleVisibility','off');
            case 4 % Line 4-5
                if options.reheat==1
                    p=linspace(p4,p5,n_pt);
                    s=linspace(s4,s5,n_pt);
                    t=zeros(size(s));
                    t(1)=t4;
                    for k=2:n_pt-1
                        t(k)=XSteam('T_ps',p(k),s(k));
                    end
                    t(end)=t5;
                    plot(s,t,'k','HandleVisibility','off');
                end
            case 5 % Line 5-6 (Expension and feed-heaters)
                n_pt=4*n_pt;
                p=linspace(p5,p6(1),n_pt);
                s=zeros(size(p));
                t=zeros(size(p));
                t(1)=t5;
                s(1)=s5;
                for k=2:n_pt-1
                    [T6,~,~,S6,~] = expansion(eta_SiT(end),...
                        t5,p5,h5,s5,x5,p(k),0);
                    t(k)=T6(end);
                    s(k)=S6(end);
                end
                t(end)=t6(1);
                s(end)=s6(1);
                plot(s,t,'k','HandleVisibility','off');
            case 6 % Line 6-7 (Heat exchangers)
                for m=1:length(t6)
                    s=linspace(s6(m),s7(m),n_pt);
                    t=zeros(size(s));
                    t(1)=t6(m);
                    for k=2:n_pt-1
                        t(k)=XSteam('t_ps',p6(m),s(k));
                    end
                    t(end)=t7(m);
                    plot(s,t,'k','HandleVisibility','off');
                end
            case 7 % Line 7-8 (Compression after the condenser and after the drum)
                p=linspace(p7(1),p8,n_pt);
                s=zeros(size(p));
                t=zeros(size(p));
                s(1)=s7(1);
                t(1)=t7(1);
                for k=2:n_pt-1
                    [t(k),~,s(k),~,~] = compression(eta_SiC,p(k),h7(1),s7(1));
                end
                t(end)=t8;
                s(end)=s8;
                plot(s,t,'k','HandleVisibility','off');
                % Pump after the drum
                if options.drumFlag==1
                    p=linspace(p7(ind_drum),p9(ind_drum),n_pt);
                    s=zeros(size(p));
                    t=zeros(size(p));
                    t(1)=t7(ind_drum);
                    s(1)=s7(ind_drum);
                    for k=2:n_pt-1
                        [t(k),~,s(k),~,~] = compression(eta_SiC,p(k),h7(ind_drum),s7(ind_drum));
                    end
                    t(end)=t9(ind_drum);
                    s(end)=s9(ind_drum);
                    plot(s,t,'k','HandleVisibility','off');
                end
            case 8 % Line 8-9o
                t=linspace(t8,t9(1),n_pt);
                s=zeros(size(t));
                s(1)=s8;
                for k=2:n_pt-1
                    s(k)=XSteam('s_pT',p8,t(k));
                end
                s(end)=s9(1);
                plot(s,t,'k','HandleVisibility','off');
                
            case 9 % Line 9o-1
                if options.drumFlag==1
                    t=linspace(t9(1),t9(ind_drum-1),n_pt);
                else
                    t=linspace(t9(1),t1,n_pt);
                end
                s=zeros(size(t));
                s(1)=s9(1);
                for k=2:n_pt
                    s(k)=XSteam('s_pT',p9(1),t(k));
                end
                plot(s,t,'k','HandleVisibility','off');
                if options.drumFlag==1 % After the drum
                    t=linspace(t9(ind_drum),t9(end),n_pt);
                    s=zeros(size(t));
                    s(1)=s9(ind_drum);
                    for k=2:n_pt
                        s(k)=XSteam('s_pT',p9(ind_drum),t(k));
                    end
                    plot(s,t,'k','HandleVisibility','off');
                end
                
        end
    end

    function plot_HS(ind)
        n_pt=100;
        switch ind
            case 1 % Line 1-2
                p=linspace(p1,p2,n_pt);
                h=zeros(size(p));
                s=zeros(size(p));
                s(1)=s1;
                h(1)=h1;
                for k=2:n_pt-1
                    [~,h(k),s(k),~,~] = compression(eta_SiC,p(k),h1,s1);
                end
                h(end)=h2;
                s(end)=s2;
                plot(s,h,'k','HandleVisibility','off');
            case 2 % Line 2-3
                p=linspace(p2,p3,n_pt);
                s=linspace(s2,s3,n_pt);
                h=zeros(size(s));
                h(1)=h2;
                for k=2:n_pt-1
                    h(k)=XSteam('h_ps',p(k),s(k));
                end
                h(end)=h3;
                plot(s,h,'k','HandleVisibility','off');
            case 3 % Line 3-4
                p=linspace(p3,p4,n_pt);
                h=zeros(size(p));
                s=zeros(size(p));
                h(1)=h3;
                s(1)=s3;
                for k=2:n_pt-1
                    [~,~,H4,S4,~] = expansion(eta_SiT(1),...
                        t3,p3,h3,s3,x3,p(k),0);
                    h(k)=H4(end);
                    s(k)=S4(end);
                end
                h(end)=h4;
                s(end)=s4;
                plot(s,h,'k','HandleVisibility','off');
            case 4 % Line 4-5
                if options.reheat==1
                    p=linspace(p4,p5,n_pt);
                    s=linspace(s4,s5,n_pt);
                    h=zeros(size(s));
                    h(1)=h4;
                    for k=2:n_pt-1
                        h(k)=XSteam('h_ps',p(k),s(k));
                    end
                    h(end)=h5;
                    plot(s,h,'k','HandleVisibility','off');
                end
            case 5 % Line 5-6 (Expension and feed-heaters)
                n_pt=4*n_pt;
                p=linspace(p5,p6(1),n_pt);
                s=zeros(size(p));
                h=zeros(size(p));
                h(1)=h5;
                s(1)=s5;
                for k=2:n_pt-1
                    [~,~,H6,S6,~] = expansion(eta_SiT(end),...
                        t5,p5,h5,s5,x5,p(k),0);
                    h(k)=H6(end);
                    s(k)=S6(end);
                end
                h(end)=h6(1);
                s(end)=s6(1);
                plot(s,h,'k','HandleVisibility','off');
            case 6 % Line 6-7 (Heat exchangers)
                for m=1:length(t6)
                    s=linspace(s6(m),s7(m),n_pt);
                    h=zeros(size(s));
                    h(1)=h6(m);
                    for k=2:n_pt-1
                        h(k)=XSteam('h_ps',p6(m),s(k));
                    end
                    h(end)=h7(m);
                    plot(s,h,'k','HandleVisibility','off');
                end
            case 7 % Line 7-8 (Compression after the condenser and after the drum)
                p=linspace(p7(1),p8,n_pt);
                s=zeros(size(p));
                h=zeros(size(p));
                s(1)=s7(1);
                h(1)=h7(1);
                for k=2:n_pt-1
                    [~,h(k),s(k),~,~] = compression(eta_SiC,p(k),h7(1),s7(1));
                end
                h(end)=h8;
                s(end)=s8;
                plot(s,h,'k','HandleVisibility','off');
                % Pump after the drum
                if options.drumFlag==1
                    p=linspace(p7(ind_drum),p9(ind_drum),n_pt);
                    s=zeros(size(p));
                    h=zeros(size(p));
                    h(1)=h7(ind_drum);
                    s(1)=s7(ind_drum);
                    for k=2:n_pt-1
                        [~,h(k),s(k),~,~] = compression(eta_SiC,p(k),h7(ind_drum),s7(ind_drum));
                    end
                    h(end)=h9(ind_drum);
                    s(end)=s9(ind_drum);
                    plot(s,h,'k','HandleVisibility','off');
                end
            case 8 % Line 8-9o
                t=linspace(t8,t9(1),n_pt);
                s=zeros(size(t));
                h=zeros(size(t));
                h(1)=h8;
                s(1)=s8;
                for k=2:n_pt-1
                    s(k)=XSteam('s_pT',p8,t(k));
                    h(k)=XSteam('h_pT',p8,t(k));
                end
                s(end)=s9(1);
                plot(s,h,'k','HandleVisibility','off');
                
            case 9 % Line 9o-1
                if options.drumFlag==1
                    t=linspace(t9(1),t9(ind_drum),n_pt);
                else
                    t=linspace(t9(1),t1,n_pt);
                end
                s=zeros(size(t));
                h=zeros(size(t));
                s(1)=s9(1);
                h(1)=h9(1);
                for k=2:n_pt
                    s(k)=XSteam('s_pT',p9(1),h(k));
                    h(k)=XSteam('h_pT',p9(1),h(k));
                end
                plot(s,h,'k','HandleVisibility','off');
                if options.drumFlag==1 % After the drum
                    t=linspace(t9(ind_drum),t9(end),n_pt);
                    h=zeros(size(t));
                    s=zeros(size(h));
                    s(1)=s9(ind_drum);
                    h(1)=h9(ind_drum);
                    for k=2:n_pt
                        s(k)=XSteam('s_pT',p9(ind_drum),h(k));
                        h(k)=XSteam('h_pT',p9(ind_drum),h(k));
                    end
                    plot(s,h,'k','HandleVisibility','off');
                end
                
                
        end
        
    end
end

