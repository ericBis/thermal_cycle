function [ETA XMASSFLOW DATEN DATEX DAT MASSFLOW COMBUSTION FIG] = ST(P_e,options,display)
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
    display = 0;
    if nargin<2
        options = struct();
        if nargin<1
            P_e = 35; % [kW] Puissance �nerg�tique de l'installation
        end
        options= struct;
        options.nsout=2; %   [-] : Number of feed-heating
        options.reheat=0; %    [-] : Number of reheating
        options.T_max=520; %     [�C] : Maximum steam temperature
        options.T_cond_out=30; %[�C] : Condenseur cold outlet temperature
        options.p3_hp=4e6; %     [bar] : Maximum pressure
        options.drumFlag=0; %   [-] : if =1 then drum if =0 => no drum.
        options.eta_mec=0.98; %    [-] : mecanic efficiency of shafts bearings
        % options.comb is a structure containing combustion data :
        options.comb.Tmax=0; %      [�C] : maximum combustion temperature
        options.comb.lambda=0; %    [-] : air excess
        options.comb.x=0; %         [-] : the ratio O_x/C. Example 0.05 in CH_1.2O_0.05
        options.comb.y=0; %         [-] : the ratio H_y/C. Example 1.2 in CH_1.2O_0.05
        options.T_exhaust=0; %  [�C] : Temperature of exhaust gas out of the chimney
        options.p_3=1050e3; %        [-] : High pressure after last reheating
        options.x4=0.91; %         [-] : Vapor ratio [gaseous/liquid] (in french : titre)
        options.T_0=15; %        [�C] : Reference temperature
        options.TpinchSub=0; %  [�C] : Temperature pinch at the subcooler
        options.TpinchEx=0; %   [�C] : Temperature pinch at a heat exchanger
        options.TpinchCond=3; % [�C] : Temperature pinch at condenser
        options.Tdrum=0; %      [�C] : minimal drum temperature
        options.eta_SiC=0.85; %     [-] : Isotrenpic efficiency for compression
        options.eta_SiT=0.88; %     [-] : Isotrenpic efficiency for Turbine. It can be a vector of 2 values :
        %             	             eta_SiT(1)=eta_SiT_HP,eta_SiT(2)=eta_SiT_others
    end
end


% Exemple of how to use (isfield' to check if an option has been given (or
% not)
if isfield(options,'T_0')
    T_0 = options.T_0+273.15;
else
    T_0 = 288.15;  % [�C]
end

%
nsout=options.nsout;
eta_SiT=options.eta_SiT;

%Prealocations
t6=zeros(nsout+1,1);
p6=zeros(nsout+1,1);
h6=zeros(nsout+1,1);
s6=zeros(nsout+1,1);
x6=zeros(nsout+1,1);
x6=zeros(nsout+1,1);

h0=XSteam('hL_T',T_0-273.15);
s0=XSteam('sL_T',T_0-273.15);

t3=options.T_max;
p3=options.p3_hp*1e-5; % [bar]
h3=XSteam('h_pT',p3,t3);
s3=XSteam('s_pT',p3,t3);
x3=nan;
e3=(h3-h0)-T_0*(s3-s0);

%Etat 4
p4=options.p_3*1e-5;
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
    x5=nan;
else
    t5=t4;
    p5=p4;
    h5=h4;
    s5=s4;
    x5=x4;
    e5=e4;
end

%Etat 7
t7=30;%options.T_cond_out+options.TpinchCond;
p7=XSteam('psat_T',t7);
h7=XSteam('hL_T',t7);
s7=XSteam('sL_T',t7);
e7=(h7-h0)-T_0*(s7-s0);
x7=0;

%Etat(s) 6
[t_i,p_i,h_i,s_i,x_i] = expansion(eta_SiT(end),t5,p5,h5,s5,x5,p7,nsout-1);
t6=flip(t_i);
p6=flip(p_i);
h6=flip(h_i);
s6=flip(s_i);
x6=flip(x_i);
e6=(h6-h0)-T_0*(s6-s0);

%Etat 8
p8=44; % [bar] ATTENTION VALEUR TEMPORERE
h8s=XSteam('h_ps',p8,s7);
h8=h7+1/options.eta_SiC*(h8s-h7);
s8=XSteam('s_ph',p8,h8);
x8=XSteam('x_ph',p8,h8);
t8=XSteam('T_ph',p8,h8);
e8=(h8-h0)-T_0*(s8-s0);

% Plots
if display ==1 
    %Cloche T-S
    Tk=374.15; %[�C] Point triple de l'eau
    T=linspace(0,Tk);
    sL_T=zeros(size(T));
    sV_T = zeros(size(T));
    hL_T=zeros(size(T));
    hV_T = zeros(size(T));
    for i=1:length(T)
        sL_T(i)=XSteam('sL_T',T(i));
        sV_T(i)=XSteam('sV_T',T(i));
        hL_T(i) = XSteam('hL_T',T(i));
        hV_T(i) = XSteam('hV_T',T(i));
    end    
    
    % 1-2
    t1_2=linespace(p1,p8,50);
    s1_2=zeros(size(t1_2));
    h1_2=zeros(size(t1_2));
    for i=1:length(t1_2)
        p2_p=p1*p3/p4;
        h2s_p=XSteam('h_ps',h2s_p,s1);
        h1_2(i)=h1+1/options.eta_SiC*(h2s_p-h1);
        s1_2(i)=XSteam('s_ph',p8,h8);
        t2=XSteam('T_ph',p8,h8);

    end

    
end

function dhdp = dComp(p,h)
    v=XSteam('v_ph',p,h);
    dhdp=v/options.SiC;
end
function dhdp= dTurb(p,h,SiT)
    v=XSteam('v_ph',p,h);
    dhdp=v*SiT;
end

% Funtions
    function [ti,pi,hi,si,xi] = expansion(eta_SiT,t_in,p_in,h_in,s_in,x_in,p_out,nsout)
        % ti,pi,hi,si,xi: state of the feed-heatings
        
        % "main" expansion
        s_outs=s_in;
        h_outs=XSteam('h_ps',p_out,s_outs);
        h_out= h_in + eta_SiT*(h_outs - h_in);
        t_out= XSteam('t_ph',p_out,h_out);
        s_out= XSteam('s_ph',p_out,h_out);
        x_out= XSteam('x_ps',p_out,s_out);
        
        % Feed-heatings
        step_h=(h_out-h_in)/(nsout+1);
        hi=h_in:step_h:h_out;
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
                [ti(i+1),pi(i+1),si(i+1),xi(i+1)]=expan_i(eta_SiT,hi(i+1),hi(1),si(1));
            end
        end
    end
              

function [t_out,p_out,s_out,x_out]=expan_i(eta_SiT,h_out,h_init,s_i)
% Compute the intermediate states of the feed-heatings

h_outs=(h_out-(1-eta_SiT)*h_init)/(eta_SiT);
s_outs=s_i;
p_out=XSteam('p_hs',h_outs,s_outs);
t_out=XSteam('t_ph',p_out,h_out);
s_out=XSteam('s_ph',p_out,h_out);
x_out = XSteam('x_ps',p_out,s_out);

if x_out < 0.88
    disp('The turbine can only work with x < 0.88')
end
end

end