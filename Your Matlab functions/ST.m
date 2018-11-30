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
        options.nsout=8; %   [-] : Number of feed-heating
        options.reheat=1; %    [-] : Number of reheating
        options.T_max=565; %     [�C] : Maximum steam temperature
        options.T_cond_out=28; %[�C] : Condenseur cold outlet temperature
        options.p3_hp=31e6; %     [bar] : Maximum pressure
        options.drumFlag=1; %   [-] : if =1 then drum if =0 => no drum.
        options.eta_mec=0.98; %    [-] : mecanic efficiency of shafts bearings
        % options.comb is a structure containing combustion data :
        options.comb.Tmax=0; %      [�C] : maximum combustion temperature
        options.comb.lambda=0; %    [-] : air excess
        options.comb.x=0; %         [-] : the ratio O_x/C. Example 0.05 in CH_1.2O_0.05
        options.comb.y=0; %         [-] : the ratio H_y/C. Example 1.2 in CH_1.2O_0.05
        options.T_exhaust=0; %  [�C] : Temperature of exhaust gas out of the chimney
        options.p_3=7e6; %        [-] : High pressure after last reheating
        options.x4=0.99; %         [-] : Vapor ratio [gaseous/liquid] (in french : titre)
        options.T_0=15; %        [�C] : Reference temperature
        options.TpinchSub=0; %  [�C] : Temperature pinch at the subcooler
        options.TpinchEx=0; %   [�C] : Temperature pinch at a heat exchanger
        options.TpinchCond=5; % [�C] : Temperature pinch at condenser
        options.Tdrum=148.7; %      [�C] : minimal drum temperature
        options.eta_SiC=0.85; %     [-] : Isotrenpic efficiency for compression
        options.eta_SiT=[0.9271 0.8874]; %     [-] : Isotrenpic efficiency for Turbine. It can be a vector of 2 values :
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

if ~isfield(options,'drumFlag')    
   options.drumFlag=1;
end

%
nsout=options.nsout;
eta_SiT=options.eta_SiT;
TpinchSub=options.TpinchSub;
TpinchEx=options.TpinchEx;
eta_SiC=options.eta_SiC;

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
p3=options.p3_hp*1e-5; % [bar]
h3=XSteam('h_pT',p3,t3);
s3=XSteam('s_pT',p3,t3);
x3=nan; % Vapeur surchauf�e
e3=(h3-h0)-T_0*(s3-s0); % [kJ/kg]

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
    p5=options.p_3*1e-5;
    h5=XSteam('h_pT',p5,t5);
    s5=XSteam('s_pT',p5,t5);
    e5=(h5-h0)-T_0*(s5-s0);
    x5=nan; % vapeur surchaf�e
else % Etat 4 = Etat 5
    t5=t4;
    p5=p4;
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
        t_sat=XSteam('Tsat_p',p6(ind_drum));
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

%Before the drum
[t_fs,h_fs,p_fs,s_fs,t9,h9,p9,s9,Xbled] = states7_9_pre_drum(h_feed,p_feed,p8,t8,h8);

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
    [t_w,h_w,p_w,s_w,t_9,h_9,p_9,s_9,Xbled2] = states7_9_post_drum(h_feed,p_feed,h7(end),s7(end));
    
    t7=[t7 ; t_w ];
    p7=[p7 ; p_w ];
    h7=[h7; h_w];
    s7=[s7 ; s_w];
    
    t9=[t9 ; t_9 ];
    p9=[p9 ; p_9 ];
    h9=[h9; h_9];
    s9=[s9 ; s_9];
    Xbled = [Xbled ; Xbled2];
end
e7= [e7 ; (h7-h0)-T_0*(s7-s0)];
e9= [e7 ; (h7-h0)-T_0*(s7-s0)];
x7=zeros(size(t7));
x9=zeros(size(t9));

%Etat 1
t1=t9(end);
p1=p9(end);
s1=s9(end);
h1=h9(end);
x1=0;

% Etat 2
p2=p3;
[t2,h2,s2,x2,e2] = compression(eta_SiC,p2,h1,s1);

% Flowrates X 

if options.drumFlag==1 
    
else
    
end



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
            disp('The turbine can only work with x < 0.88')
        end
    end

    function [t_w,h_w,p_w,s_w,t_9,h_9,p_9,s_9,Xbled] = states7_9_pre_drum(h6,p6,p8,t8,h8)
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
        
        % State sc : Between the subcooler and the Condensor
        tsc=t8+TpinchSub;
        psc=p_w(1);
        hsc=XSteam('h_pt',psc,tsc);
        
        % Computations of the X
        % They can be compbuted by resolving a linear system of type Ax=B
        % Here x= [h9_0 X_I X_II ...]
        
        A=zeros(ns+1,ns+1);
        A(end-1,1)=1;
        A(end,1)=-1;
        A(end,2:end)=hsc-h_w(1);
        
        %Filling A
        k=ns+1;
        for i=2:ns
            A(i,k:end)=h_w(k-1)-h_w(k-2);
            k=k-1;
        end
        
        k=ns;
        for i=1:ns
            A(i,k+1)=h6(k)-h_w(k);
            k=k-1;
        end
        
        B=zeros(ns+1,1);
        B(end-1)=h_9(2);
        B(end)=-h8;
        
        k=ns;
        for i=1:ns-1
            B(i)=h_9(k+1)-h_9(k);
            k=k-1;
        end
        
        x=A\B;
        h_9(1)=x(1);
        t_9(1)=XSteam('t_ph',p8,h_9(1));
        s_9(1)=XSteam('s_ph',p8,h_9(1));
        Xbled=x(2:end);
    end

function [t_w,h_w,p_w,s_w,t_9,h_9,p_9,s_9,Xbled] = states7_9_post_drum(h6,p6,h_bef,s_bef)
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
        
        for i=2:ns
            h_9(i)=XSteam('h_pT',p_9(i),t_9(i));
            s_9(i)=XSteam('s_pT',p_9(i),t_9(i));
        end
        
        % Compression at Pb
        [t_9(1),h_9(1),s_9(1),~,~] = compression(options.eta_SiC,p_9(1),h_bef,s_bef);
        

        % Computations of the X
        % They can be computed by resolving a linear system of type Ax=B
        % Here x= [ X_i X_i+1 ... X_N]
        
        A=zeros(ns,ns);
        
        %Filling A
        k=ns;
        for i=2:ns
            A(i,k:end)=h_w(k)-h_w(k-1);
            k=k-1;
        end
        
        k=ns;
        for i=1:ns
            A(i,k)=h6(k)-h_w(k);
            k=k-1;
        end
        
        B=zeros(ns,1);         
        k=ns;
        
        for i=1:ns
            B(i)=h_9(k+1)-h_9(k);
            k=k-1;
        end
        
        Xbled=A\B;
        t_9(1)=XSteam('t_ph',p_9(1),h_9(1));
        s_9(1)=XSteam('s_ph',p_9(1),h_9(1));
end

    function [X3, X5, X6, Xd] = X_flowrates(h6,h7,h9,Xbled1,Xbled2,drumFlag,ind_drum)
        %Computes the X_flowrates 
        % Xbled1 : X_i before the drum
        % Xbled2: X_i after the drum
        
        if drumFlag == 1
            A=zeros(4,4);           
            A(1,1:end)=[sum(Xbled2)*(h7(ind_drum+1)-h7(ind_drum))...
                0 0 h6(ind_drum)-h7(ind_drum)];
            A(2,1:end)=[sum(Xbled2)-Xbled2(end) -1 0 1];
            A(3,1:end)=[0 0 1 0];
            A(4,1:end)=[1-Xbled2(end) -1 0 0];
            
            B = [h7(ind_drum) - h9(ind_drum-1) ; -1 ; 1-sum(Xbled1) ; 0];
            
            x=A\B;            
            X3=x(1);
            X5=x(2);
            X6=x(3);
            Xd=x(4); % X_flowrate in the drum
        else
            A=zeros(3,3);            
            A(1,1:end)=[1 -1 0];
            A(2,1:end)=[0 1 -1];
            A(3,1:end)=[0 0 1];
            
            B=[Xbled1(end) ; sum(Xbled1)-Xbled1(end) ; 1 - sum(Xbled1)];
            
            x=A\B;            
            X3=x(1);
            X5=x(2);
            X6=x(3);
        end
        
        x
        
    end

    function [LHV,h_f,e_f,lambda,cp_g,comp_fum,h_fum,h_air,ma1] = combustion(comb)
        x=comb.x;
        y=comb.y;
        ma1=(1+(y-2*x)/4)*(32+3.76*28.15)/(12.01+1.008*y+16*x);
        LHV=getLHV(y,x);
        
    end

    function LHV =getLHV(y,x)
        LHV=393400+102250*y-(x/(1+0.5*y))*(111000+102250*y)/(12+y+x*16);
    end

end

