function [ETA MASSFLOW FIG] = CCGT3P(P_eg,options,display)
% CCGT3P is a Combine cycle Gas Turbine with 2 pressure level
% CCGT3P(P_e,options,display) compute the thermodynamics states for a CCGT
% with 3 pressure level (cfr p166 english reference book) including
% combustion, exchanger and cycles. This is done based on several inputs
% (given in OPTION) and based on a given electricity production P_e.
% It returns the main results. It can as well plots graphs if input
% argument DISPLAY = true (<=> DISPLAY=1)
%
% INPUTS (some inputs can be dependent on others => only one of these 2 can
%         be activated)
% P_EG = electrical power output target for gas turbine [kW]
% OPTIONS is a structure containing :
%   -options.T0       [�C] : Reference temperature
%   -options.T_ext    [�C] : External temperature
%   -options.T_STmax  [�C] : maximum temperature on ST cycle
%   -options.eta_mec  [-] : mecanic efficiency of shafts bearings
%   -options.pdrum   [bar]: Drum pressure
%   -options.pmid    [bar]: Intermediary pressure level
%   -options.x7       [-] : Vapor ratio [gaseous/liquid] (titre)
%   -option.eta_SiC   [-] : Isotrenpic efficiency for compression
%   -option.eta_SiT   [-] : Isotrenpic efficiency for compression
%   -options.GT    [struct] : options for Gas turbine (see GT function)
% DISPLAY = 1 or 0. If 1, then the code should plot graphics. If 0, then
%          do not plot.
%
%OUPUTS :
% ETA is a vector with :
%   -eta(1)  : eta_STcyclen, cycle energy efficiency
%   -eta(2)  : eta_GTcyclen, cycle energy efficiency
%   -eta(3)  : eta_toten, overall energy efficiency
%   -eta(4)  : eta_STcyclex, cycle exegy efficiency
%   -eta(5)  : eta_GTcyclex, cycle exegy efficiency
%   -eta(6)  : eta_totex, overall exergie efficiency
%   -eta(7)  : eta_gen, Steam generator energy efficiency
%   -eta(8)  : eta_gex, Steam generator exergy efficiency
%   -eta(9)  : eta_combex, Combustion exergy efficiency
%   -eta(10) : eta_chemex, Chimney exergy efficiency (losses)
%   -eta(11) : eta_transex, Heat exchanger overall exergy efficiency
%   FYI : eta(i) \in [0;1] [-]
% MASSFLOW is a vector containing :
%   -massflow(1) [kg/s]: water massflow at high pressure turbine inlet
%   -massflow(2) [kg/s]: water massflow at medium pressure turbine inlet
%   -massflow(3) [kg/s]: water massflow at low pressure turbine inlet
%   -massflow(4) [kg/s]: air massflow at gas turbine inlet
%   -massflow(5) [kg/s]: combustible massflow
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

%% Your work
if nargin<3
    display=1;
    if nargin<2
        options=struct();
        if nargin<1
            P_eg=225e3;%100MW
        end
        %%%%%DESCRIPTION DE LA STRUCTURE%%%%%
%           P_eg = 225e3; %[kW]
%           options=struct;
%           options.T0 =15; %[�C]:Reference temperature
%           options.T_ext=15; %[�C]:External temperature
%           options.T_STmax=565; %[�C]:maximum temperature on ST cycle
%           options.eta_mec=0.96; %[-]:mecanic efficiency of shafts bearings
%           options.pdrum=4; %[bar]:Drum pressure
%           options.pmid=28; %[bar]:Intermediary pressure level
%           options.x7=0.95;
%           options.eta_SiC=0.857;
%           options.eta_SiT=0.857;
%           options.GT=struct;%[struct]
%                    options.GT.k_mec=0.015;
%                    options.GT.T_0=15;
%                    options.GT.T_ext=15;
%                    options.GT.r=15;
%                    options.GT.k_cc=0.95;
%                    options.GT.T_3=1250;
%                    options.GT.eta_PiC=0.9;
%                    options.GT.eta_PiT=0.9;
%                    options.GT.NTU=0;
%                    options.GT.ER=0;
    end
end

if isfield(options,'T_0')
    T_0 = options.T_0;
else
    T_0 = 15; % [�C]
end

if ~isfield(options,'T0')
    options.T0=15;
end

if ~isfield(options,'T_ext')
    options.T_ext=15;
end

if ~isfield(options,'T_STmax')
    options.T_STmax=565;
end

if ~isfield(options,'eta_mec')
    options.eta_mec=0.96;
end

if ~isfield(options,'pdrum')
    options.pdrum=4;
end

if ~isfield(options,'pmid')
    options.pmid=28;
end

if ~isfield(options,'x7')
    options.x7=0.95;
end

if ~isfield(options,'eta_SiT')
    options.eta_SiT=0.857;
end

if ~isfield(options,'eta_SiC')
    options.eta_SiC=0.857;
end

if ~isfield(options,'GT')
    options.GT=struct;
end

if ~isfield(options.GT,'k_mec')
    options.GT.k_mec=0.015;
end

if ~isfield(options.GT,'T_0')
    options.GT.T_0=15;
end

if ~isfield(options.GT,'T_ext')
    options.GT.T_ext=15;
end

if ~isfield(options.GT,'r')
    options.GT.r=15;
end

if ~isfield(options.GT,'k_cc')
    options.GT.k_cc=0.95;
end

if ~isfield(options.GT,'T_3')
    options.GT.T_3=1250;
end

if ~isfield(options.GT,'eta_PiC')
    options.GT.eta_PiC=0.9;
end

if ~isfield(options.GT,'eta_PiT')
    options.GT.eta_PiT=0.9;
end

if ~isfield(options.GT,'NTU')
    options.GT.NTU=0;
end

if ~isfield(options.GT,'ER')
    options.GT.ER=0;
end

if ~isfield(options,'GraphCompute')
    options.GraphCompute=1;
end

%%%%OUTPUT GT%%%%
[ETA, DATEN, DATEX, DAT, MASSFLOW, COMBUSTION,FIG] = GT(P_eg,options.GT,display);
ETAGT=ETA;
MASSFLOWGT=MASSFLOW;

%%%%OUTPUT VARIABLE CCGT%%%%
ETA=ones(1,11);
MASSFLOW=ones(1,5);
MASSFLOW(4)=MASSFLOWGT(1);%MFAir
MASSFLOW(5)=MASSFLOWGT(2);%MFComb
MFumee=MASSFLOWGT(3);%MFumee

%%%%FUMEE EN 4G%%%%
T_4g=DAT(1,4);%[�C]
h_4g=DAT(3,4);
E_4g=DAT(5,4);
CpMoyen_4g=h_4g/T_4g;%Dans la ccgt on considere que le cpmoyen est constant

%%%%Calcul des etats de la CCGT%%%%
%Variables Etat 0:
H0=XSteam('hL_T',options.T0);
T0=273.15+options.T0;%[K]
S0=XSteam('sL_T',options.T0);
H_0g=DAT(3,1);
S_0g=DAT(4,1);

%Pinch et DeltaEchangeur
DeltaTa=T_4g-options.T_STmax;%[�C] difference de temperature entre T4g et T3
DeltaPinch=8;%[�C] difference de temperature entre fumees et eau co_courante
DeltaCondenseur=21;%[�C] difference de temperature entre riviere et vapeur

%%%Etats du cycle
T3=options.T_STmax;%[�C] difference fumee en sortie de GT et eau

%L'etat 5 est entierement connu sur base des inputs
T5=T3;%[�C] IMPOSEE
P5=options.pmid;%[bar]
H5=XSteam('h_pT',P5,T5);
S5=XSteam('s_pT',P5,T5);
E5=(H5-H0)-T0*(S5-S0);

%L'etat 9 est determine � partir de l'etat 5
P9_prime=P5;
P9_seconde=P5;
T9_prime=XSteam('Tsat_p',P9_prime);
T9_seconde=T9_prime;
X9_prime=0;
X9_seconde=1;
H9_prime=XSteam('hL_p',P9_prime);
H9_seconde=XSteam('hV_p',P9_seconde);
S9_prime=XSteam('sL_p',P9_prime);
S9_seconde=XSteam('sV_p',P9_seconde);
E9_prime=(H9_prime-H0)-T0*(S9_prime-S0) ;
E9_seconde=(H9_seconde-H0)-T0*(S9_seconde-S0) ;

%L'etat 8 est determinee a partir de l'etat 6
P6=options.pdrum;%[bar] IMPOSEE
P8=P6;
T8=T9_prime;
H8=XSteam('h_pT',P8,T8);
S8=XSteam('s_pT',P8,T8);
E8=(H8-H0)-T0*(S8-S0);
P8_prime=P6;
P8_seconde=P6;
T8_prime=XSteam('Tsat_p',P8_prime);
T8_seconde=T8_prime;
X8_prime=0;
X8_seconde=1;
H8_prime=XSteam('hL_p',P8_prime);
H8_seconde=XSteam('hV_p',P8_seconde);
S8_prime=XSteam('sL_p',P8_prime);
S8_seconde=XSteam('sV_p',P8_seconde);
E8_prime=(H8_prime-H0)-T0*(S8_prime-S0) ;
E8_seconde=(H8_seconde-H0)-T0*(S8_seconde-S0) ;

%L'etat 6 ne depend que de l'etat 5
S6s=S5;
P6s=P6;
H6s=XSteam('h_ps',P6s,S6s);
H6=H5-options.eta_SiT*(H5-H6s);
S6=XSteam('s_ph',P6,H6);
T6=XSteam('T_hs',H6,S6);
E6=(H6-H0)-T0*(S6-S0);

%Definir ici l'etat 10 et 3
T10_prime=T6;
T10_seconde=T10_prime;
X10_prime=0;
X10_seconde=1;
P10_prime=XSteam('psat_T',T10_prime);%[bar]
P10_seconde=P10_prime;%[bar]
H10_prime=XSteam('hL_p',P10_prime);
H10_seconde=XSteam('hV_p',P10_seconde);
S10_prime=XSteam('sL_p',P10_prime);
S10_seconde=XSteam('sV_p',P10_seconde);
E10_prime=(H10_prime-H0)-T0*(S10_prime-S0) ;
E10_seconde=(H10_seconde-H0)-T0*(S10_seconde-S0) ;

%L'etat 3 est defini a partir de 10
P3=P10_seconde;
H3=XSteam('h_pT',P3,T3);%[kJ/kg]
S3=XSteam('s_pT',P3,T3);%[kJ/kgK]
E3=(H3-H0)-T0*(S3-S0);%[kJ/kg]

P4=P5;
P4s=P4;
S4s=S3;
H4s=XSteam('h_ps',P4s,S4s);
H4=H3-options.eta_SiT*(H3-H4s);
S4=XSteam('s_ph',P4,H4);
T4=XSteam('T_ps',P4,S4);
E4=(H4-H0)-T0*(S4-S0);

%L'etat 9 est determiner a partur de 5 et 6
T9=T6;
P9=P5;
H9=XSteam('h_pT',P9,T9);
S9=XSteam('s_pT',P9,T9);
E9=(H9-H0)-T0*(S9-S0);

X7=options.x7;
T7=DeltaCondenseur+options.T0;
P7=XSteam('psat_T',T7);
H7=XSteam('h_Tx',T7,X7);
S7=XSteam('s_ph',P7,H7);
E7=(H7-H0)-T0*(S7-S0);

X1=0;
T1=T7;
P1=P7;
H1=XSteam('hL_T',T1);
S1=XSteam('sL_T',T1);
E1=(H1-H0)-T0*(S1-S0);

P2=P6;
P2s=P2;
S2s=S1;
H2s=XSteam('h_ps',P2s,S2s);
H2=H1+(1/options.eta_SiC)*(H2s-H1);
S2=XSteam('s_ph',P2,H2);
T2=XSteam('T_ps',P2,S2);
E2=(H2-H0)-T0*(S2-S0);

Temperatures=[T1;T2;T3;T4;T5;T6;T7;T8;T8_prime;T8_seconde;T9;T9_prime;T9_seconde;T10_prime;T10_seconde];
Pressions=[P1;P2;P3;P4;P5;P6;P7;P8;P8_prime;P8_seconde;P9;P9_prime;P9_seconde;P10_prime;P10_seconde];
Enthalpies=[H1;H2;H3;H4;H5;H6;H7;H8;H8_prime;H8_seconde;H9;H9_prime;H9_seconde;H10_prime;H10_seconde];
Entropies=[S1;S2;S3;S4;S5;S6;S7;S8;S8_prime;S8_seconde;S9;S9_prime;S9_seconde;S10_prime;S10_seconde];
Exergie=[E1;E2;E3;E4;E5;E6;E7;E8;E8_prime;E8_seconde;E9;E9_prime;E9_seconde;E10_prime;E10_seconde];
MatriceEtats=ones(length(Temperatures),5);
MatriceEtats(:,1)=Temperatures;
MatriceEtats(:,2)=Pressions;
MatriceEtats(:,3)=Enthalpies;
MatriceEtats(:,4)=Entropies;
MatriceEtats(:,5)=Exergie;

%%%%DONNEES UTILE DE LA GT POUR LA DETERMINATION DES MASSFLOW%%%%
h_HPg=CpMoyen_4g*(T10_prime+DeltaPinch);
h_IPg=CpMoyen_4g*(T9_prime+DeltaPinch);
h_LPg=CpMoyen_4g*(T8_prime+DeltaPinch);

%%%%MASSFLOW de la CCGT%%%%
massflowWater_HP=1;
massflowWater_IP=1;
massflowWater_LP=1;
x0=[massflowWater_HP,massflowWater_IP,massflowWater_LP];
funFlow=@CalculationFlow;
mHPmMIPmLP=fsolve(@(x)funFlow(x,MFumee,h_4g,h_HPg,h_IPg,h_LPg,Enthalpies),x0);
MASSFLOW(1)=mHPmMIPmLP(1);%Complete le vecteur avec les bons massflow
MASSFLOW(2)=mHPmMIPmLP(2);%Complete le vecteur avec les bons massflow
MASSFLOW(3)=mHPmMIPmLP(3);%Complete le vecteur avec les bons massflow


%%%%SORTIE DES GAZ DE LA CHEMINEE%%%%
h_5g=h_LPg-((sum(MASSFLOW(1:3)))/MFumee)*(Enthalpies(9)-Enthalpies(2));
T_5g=h_5g/CpMoyen_4g;%[�C]
S_5g=CpMoyen_4g*log((T_5g+273.15)/273.15);
E_5g=(h_5g-H_0g)-T0*(S_5g-S_0g);

%%%%PERTE PUISSANCE EXERGETIQUE en kW%%%%
IrrevExergCombustion=(MASSFLOW(5)*COMBUSTION.e_c+MASSFLOW(4)*DAT(5,2))-(MFumee)*DAT(5,3);
PuissanceEffectiveGT=P_eg;
PuissanceExergPrimaire=(MASSFLOW(5)*COMBUSTION.e_c);

Q_point=MFumee*(h_4g-h_5g);
TIm=(T_5g-T_4g)/log((T_5g+273.15)/(T_4g+273.15));%[K]
TIIm=(T3-T2)/log((T3+273.15)/(T2+273.15));%[K]
IrrevExergTransfThermique=Q_point*T0*(TIm-TIIm)/(TIm*TIIm);

PerteExergCheminee=MFumee*E_5g;
PerteExergCondenseur=sum(MASSFLOW(1:3))*(E7-E1);

irrevCompresseurGT=T0*MASSFLOW(4)*(DAT(4,2)-DAT(4,1));%AirGT
irrevTurbineGT=T0*MFumee*(DAT(4,4)-DAT(4,3));%FumeeGT
irrevTurbineST_HP=T0*MASSFLOW(1)*(S4-S3);%mvHP dans la turbine HP
irrevTurbineST_IP=T0*sum(MASSFLOW(1:2))*(S6-S5);%mvHP+mvIP dans la turbine IP
irrevTurbineST_LP=T0*sum(MASSFLOW(1:3))*(S7-S6);%mvHP+mvIP+mvLP dans la turbine LP
IrrevComplexeRotoriqueCCGT=irrevCompresseurGT+irrevTurbineGT+irrevTurbineST_HP+irrevTurbineST_IP+irrevTurbineST_LP;

PmST_turbines=MASSFLOW(1)*(H3-H4)+sum(MASSFLOW(1:2))*(H5-H6)+sum(MASSFLOW(1:3))*(H6-H7);
PmST=PmST_turbines;
PuissanceEffectiveST=PmST*options.eta_mec;
PerteMecaST=PmST-PuissanceEffectiveST;
PerteMecaGT=DATEN(1);
PerteMecaCCGT=PerteMecaST+PerteMecaGT;

%%%%PERTE PUISSANCE ENERGETIQUE en kW%%%%
PerteEnergCheminee=MFumee*(h_5g-DAT(3,1));
PerteEnergCondenseur=sum(MASSFLOW(1:3))*(H7-H1);
PuissanceEnergPrimaire=(MASSFLOW(5)*COMBUSTION.LHV);


%Plots
if options.GraphCompute==1
    if display ==1
        disp='on';
    else
        disp='off';
    end
    
    
%     %Cloche T-S
%     S_plot=linspace(0,9.156,300);
%     T_HS=zeros(size(S_plot));
%     for c=1:length(T_HS)
%         T_HS(c)=XSteam('Tsat_s',S_plot(c));
%     end
%     
%     H_HS=zeros(size(S_plot));
%     [~,ind_tmax]=max(T_HS);
%     for c=1:ind_tmax
%         H_HS(c)=XSteam('hL_T',T_HS(c));
%     end
%     for c=ind_tmax+1:length(T_HS)
%         H_HS(c)=XSteam('hV_T',T_HS(c));
%     end
        
    FIG = gobjects(4,1);
    
%     FIG(1)=figure('visible',disp);
%     hold on;
%     plot(S_plot,T_HS);
%     title('CCGT T-s diagram');
%     ylabel('Temperature [�C]');
%     xlabel('Entropy [kJ/kg/K]');
%     for i=1:17
%         plot_TS(i);
%     end
%     
%     FIG(2)=figure('visible',disp);
%     hold on;
%     plot(S_plot,H_HS);
%     title('CCGT h-s diagram');
%     ylabel('Enthalpy [kJ/kg]');
%     xlabel('Entropy [kJ/kg/K]');
%     for i=1:17
%         plot_HS(i);
%     end
    
%     %%%%PIE CHART EXERG%%%%
%     FIG(3)=figure('visible',disp);
%     dataExerg = [IrrevExergCombustion,IrrevExergTransfThermique,PerteExergCheminee,IrrevComplexeRotoriqueCCGT,PerteExergCondenseur,PerteMecaCCGT,PuissanceEffectiveST,PuissanceEffectiveGT]*(10^(-3));
%     labelsExerg = {strcat('Irrev. comb.:',num2str(dataExerg(1)),'MW'),strcat('Irrev. transf.:',num2str(dataExerg(2)),'MW'),strcat('Perte cheminee:',num2str(dataExerg(3)),'MW'),strcat('Irrev. rotor CCGT:',num2str(dataExerg(4)),'MW'),strcat('Perte condenseur:',num2str(dataExerg(5)),'MW'),strcat('Perte meca CCGT:',num2str(dataExerg(6)),'MW'),strcat('Puiss. TV:',num2str(dataExerg(7)),'MW'),strcat('Puiss. TG:',num2str(dataExerg(8)),'MW')};
%     pie(dataExerg,labelsExerg);
%     title(strcat('Puiss. exerg�tique primaire:',num2str(PuissanceExergPrimaire*(10^(-3))),'MW'));
%     
%     %%%%PIE CHART ENERG%%%%
%     FIG(4)=figure('visible',disp);
%     dataEnerg = [PerteEnergCheminee,PerteEnergCondenseur,PerteMecaCCGT,PuissanceEffectiveST,PuissanceEffectiveGT]*(10^(-3));
%     labelsEnerg = {strcat('Perte cheminee:',num2str(dataEnerg(1)),'MW'),strcat('Perte condenseur:',num2str(dataEnerg(2)),'MW'),strcat('Perte meca CCGT:',num2str(dataEnerg(3)),'MW'),strcat('Puiss. TV:',num2str(dataEnerg(4)),'MW'),strcat('Puiss. TG:',num2str(dataEnerg(5)),'MW')};
%     pie(dataEnerg,labelsEnerg);
%     title(strcat('Puiss. energ�tique primaire:',num2str(PuissanceEnergPrimaire*(10^(-3))),'MW'));
    
end

%%%%RENDEMENT DE LA ST, GT ET CCGT%%%%
ETA(1)=PmST/Q_point;%eta_STcyclen, cycle energy efficiency
ETA(2)=ETAGT(1);%eta_GTcyclen
ETA(3)=(PuissanceEffectiveGT+PuissanceEffectiveST)/(MASSFLOW(5)*COMBUSTION.LHV);%eta_toten
ETA(4)=PmST/(MFumee*(E_4g-E_5g));%eta_STcyclex, cycle exegy efficiency
ETA(5)=ETAGT(3);%eta_GTcyclex, cycle exegy efficiency
ETA(6)=(PuissanceEffectiveGT+PuissanceEffectiveST)/(MASSFLOW(5)*COMBUSTION.e_c);%eta_totex, overall exergie efficiency
ETA(7)=(MFumee*(h_4g-h_5g))/(MASSFLOW(5)*COMBUSTION.LHV) ;%eta_gen, Steam generator energy efficiency
ETA(9)=ETAGT(6);%eta_combex, Combustion exergy efficiency
ETA(10)=(E_4g-E_5g)/(E_4g) ;%eta_chemex, Chimney exergy efficiency (losses)
NumerateurEtaTransexSG=MASSFLOW(1)*(E3)-(sum(MASSFLOW(1:3)))*E2;
DenominateurEtaTransexSG=MFumee*(E_4g-E_5g);
ETA(11)=NumerateurEtaTransexSG/DenominateurEtaTransexSG;%eta_transex, Heat exchanger overall exergy efficiency
ETA(8)=ETA(10)*ETA(11)* ETA(9); %eta_gex, Steam generator exergy efficiency

    function systemMassFlow = CalculationFlow(MassflowIterInitiale,MFumee,h_4g,h_HPg,h_IPg,h_LPg,Enthalpies)
        systemMassFlow=ones(3,1);
        systemMassFlow(1,1)= -MFumee*(h_4g-h_HPg)+MassflowIterInitiale(1)*(Enthalpies(3)-Enthalpies(14)+Enthalpies(5)-Enthalpies(4))+MassflowIterInitiale(2)*(Enthalpies(5)-Enthalpies(11));
        systemMassFlow(2,1)= -MFumee*(h_HPg-h_IPg)+MassflowIterInitiale(1)*(Enthalpies(14)-Enthalpies(12))+MassflowIterInitiale(2)*(Enthalpies(11)-Enthalpies(12))+MassflowIterInitiale(3)*(Enthalpies(6)-Enthalpies(8));
        systemMassFlow(3,1)= -MFumee*(h_IPg-h_LPg)+MassflowIterInitiale(1)*(Enthalpies(12)-Enthalpies(9))+MassflowIterInitiale(2)*(Enthalpies(12)-Enthalpies(9))+MassflowIterInitiale(3)*(Enthalpies(8)-Enthalpies(9));
    end

    function [t2,h2,s2,x2,e2] = compression(eta_SiC,p2,h1,s1)
        % Adiabatic compression
        h_2s=XSteam('h_ps',p2,s1);
        h2=h1+1/eta_SiC*(h_2s-h1);
        s2=XSteam('s_ph',p2,h2);
        x2=XSteam('x_ph',p2,h2);
        t2=XSteam('T_ph',p2,h2);
        e2=(h2-H0)-T0*(s2-S0);
    end

    function [t_out,p_out,s_out,x_out]=expan_inter(eta_SiT,h_out,h_init,s_i)
        
        h_outs=(h_out-(1-eta_SiT)*h_init)/(eta_SiT); % ISentropic expansion
        s_outs=s_i;
        p_out=XSteam('p_hs',h_outs,s_outs);
        t_out=XSteam('t_ph',p_out,h_out);
        s_out=XSteam('s_ph',p_out,h_out);
        x_out = XSteam('x_ps',p_out,s_out);
        
        if x_out < 0.88
            error('The turbine can only work with x < 0.88')
        end
    end

    function [t2,h2,s2,x2,e2] = detente(eta_SiT,p2,h1,s1)
        % Adiabatic compression
        h_2s=XSteam('h_ps',p2,s1);
        h2= h1 + eta_SiT*(h_2s - h1);
        s2=XSteam('s_ph',p2,h2);
        x2=XSteam('x_ph',p2,h2);
        t2=XSteam('T_ph',p2,h2);
        e2=(h2-H0)-T0*(s2-S0);
    end

%     function plot_TS(ind)
%         n_pt=100;
%         switch ind
%             case 1 % Line 1-2
%                 p=linspace(P1,P2,n_pt);
%                 t=zeros(size(p));
%                 s=zeros(size(p));
%                 s(1)=S1;
%                 t(1)=T1;
%                 for k=2:n_pt-1
%                     [t(k),~,s(k),~,~] = compression(options.eta_SiC,p(k),H1,S1);
%                 end
%                 t(end)=T2;
%                 s(end)=S2;
%                 plot(s,t,'k','HandleVisibility','off');
%             case 2 % Line 3-4
%                 n_pt=4*n_pt;
%                 p=linspace(P3,P4,n_pt);
%                 s=zeros(size(p));
%                 t=zeros(size(p));
%                 t(1)=T3;
%                 s(1)=S3;
%                 for k=2:n_pt-1
%                     [t(k),~,s(k),~] = expan_inter(options.eta_SiT,H4,H3,S3);
%                 end
%                 t(end)=T4;
%                 s(end)=S4;
%                 plot(s,t,'k','HandleVisibility','off');
%             case 3 % Line 5-6
%                 n_pt=4*n_pt;
%                 p=linspace(P5,P6,n_pt);
%                 s=zeros(size(p));
%                 t=zeros(size(p));
%                 t(1)=T5;
%                 s(1)=S5;
%                 for k=2:n_pt-1
%                     [t(k),~,s(k),~] = expan_inter(options.eta_SiT,H6,H5,S5);
%                 end
%                 t(end)=T6;
%                 s(end)=S6;
%                 plot(s,t,'k','HandleVisibility','off');
%             case 4 % Line 6-7
%                 n_pt=4*n_pt;
%                 p=linspace(P6,P7,n_pt);
%                 s=zeros(size(p));
%                 t=zeros(size(p));
%                 t(1)=T6;
%                 s(1)=S6;
%                 for k=2:n_pt-1
%                     [t(k),~,s(k),~] = expan_inter(options.eta_SiT,H7,H6,S6);
%                 end
%                 t(end)=T7;
%                 s(end)=S7;
%                 plot(s,t,'k','HandleVisibility','off');
%             case 5 % Line 8'-9'
%                 p=linspace(P8_prime,P9_prime,n_pt);
%                 t=zeros(size(p));
%                 s=zeros(size(p));
%                 s(1)=S8_prime;
%                 t(1)=T8_prime;
%                 for k=2:n_pt-1
%                     [t(k),~,s(k),~,~] = compression(options.eta_SiC,p(k),H8_prime,S8_prime);
%                 end
%                 t(end)=T9_prime;
%                 s(end)=S9_prime;
%                 plot(s,t,'k','HandleVisibility','off');
%             case 6 % Line 9'-10'
%                 p=linspace(P9_prime,P10_prime,n_pt);
%                 t=zeros(size(p));
%                 s=zeros(size(p));
%                 s(1)=S9_prime;
%                 t(1)=T9_prime;
%                 for k=2:n_pt-1
%                     [t(k),~,s(k),~,~] = compression(options.eta_SiC,p(k),H9_prime,S9_prime);
%                 end
%                 t(end)=T10_prime;
%                 s(end)=S10_prime;
%                 plot(s,t,'k','HandleVisibility','off');
%             case 7 % Line 8'-8''
%                 p=linspace(P8_prime,P8_seconde,n_pt);
%                 s=linspace(S8_prime,S8_seconde,n_pt);
%                 t=zeros(size(s));
%                 t(1)=T8_prime;
%                 for k=2:n_pt-1
%                     t(k)=XSteam('T_ps',p(k),s(k));
%                 end
%                 t(end)=T8_seconde;
%                 plot(s,t,'k','HandleVisibility','off');
%             case 8 % Line 9'-9''
%                 p=linspace(P9_prime,P9_seconde,n_pt);
%                 s=linspace(S9_prime,S9_seconde,n_pt);
%                 t=zeros(size(s));
%                 t(1)=T9_prime;
%                 for k=2:n_pt-1
%                     t(k)=XSteam('T_ps',p(k),s(k));
%                 end
%                 t(end)=T9_seconde;
%                 plot(s,t,'k','HandleVisibility','off');
%             case 9 % Line 10'-10''
%                 p=linspace(P10_prime,P10_seconde,n_pt);
%                 s=linspace(S10_prime,S10_seconde,n_pt);
%                 t=zeros(size(s));
%                 t(1)=T10_prime;
%                 for k=2:n_pt-1
%                     t(k)=XSteam('T_ps',p(k),s(k));
%                 end
%                 t(end)=T10_seconde;
%                 plot(s,t,'k','HandleVisibility','off');
%             case 10 % Line 8''-8
%                 p=linspace(P8_seconde,P8,n_pt);
%                 s=linspace(S8_seconde,S8,n_pt);
%                 t=zeros(size(s));
%                 t(1)=T8_seconde;
%                 for k=2:n_pt-1
%                     t(k)=XSteam('T_ps',p(k),s(k));
%                 end
%                 t(end)=T8;
%                 plot(s,t,'k','HandleVisibility','off');
%             case 11 % Line 8-6
%                 p=linspace(P8,P6,n_pt);
%                 s=linspace(S8,S6,n_pt);
%                 t=zeros(size(s));
%                 t(1)=T8;
%                 for k=2:n_pt-1
%                     t(k)=XSteam('T_ps',p(k),s(k));
%                 end
%                 t(end)=T6;
%                 plot(s,t,'k','HandleVisibility','off');
%             case 12 % Line 9''-9
%                 p=linspace(P9_seconde,P9,n_pt);
%                 s=linspace(S9_seconde,S9,n_pt);
%                 t=zeros(size(s));
%                 t(1)=T9_seconde;
%                 for k=2:n_pt-1
%                     t(k)=XSteam('T_ps',p(k),s(k));
%                 end
%                 t(end)=T9;
%                 plot(s,t,'k','HandleVisibility','off');
%             case 13 % Line 10''-3
%                 p=linspace(P10_seconde,P3,n_pt);
%                 s=linspace(S10_seconde,S3,n_pt);
%                 t=zeros(size(s));
%                 t(1)=T10_seconde;
%                 for k=2:n_pt-1
%                     t(k)=XSteam('T_ps',p(k),s(k));
%                 end
%                 t(end)=T3;
%                 plot(s,t,'k','HandleVisibility','off');
%             case 14 % Line 2-8'
%                 p=linspace(P2,P8_prime,n_pt);
%                 s=linspace(S2,S8_prime,n_pt);
%                 t=zeros(size(s));
%                 t(1)=T2;
%                 for k=2:n_pt-1
%                     t(k)=XSteam('T_ps',p(k),s(k));
%                 end
%                 t(end)=T8_prime;
%                 plot(s,t,'k','HandleVisibility','off');
%             case 15 % Line 9-4
%                 p=linspace(P9,P4,n_pt);
%                 s=linspace(S9,S4,n_pt);
%                 t=zeros(size(s));
%                 t(1)=T9;
%                 for k=2:n_pt-1
%                     t(k)=XSteam('T_ps',p(k),s(k));
%                 end
%                 t(end)=T4;
%                 plot(s,t,'k','HandleVisibility','off');
%             case 16 % Line 2-7
%                 p=linspace(P1,P7,n_pt);
%                 s=linspace(S1,S7,n_pt);
%                 t=zeros(size(s));
%                 t(1)=T1;
%                 for k=2:n_pt-1
%                     t(k)=XSteam('T_ps',p(k),s(k));
%                 end
%                 t(end)=T7;
%                 plot(s,t,'k','HandleVisibility','off');
%             case 17 % Line 4-5
%                 p=linspace(P4,P5,n_pt);
%                 s=linspace(S4,S5,n_pt);
%                 t=zeros(size(s));
%                 t(1)=T4;
%                 for k=2:n_pt-1
%                     t(k)=XSteam('T_ps',p(k),s(k));
%                 end
%                 t(end)=T5;
%                 plot(s,t,'k','HandleVisibility','off');
%         end
%     end
% 
%     function plot_HS(ind)
%         n_pt=100;
%         switch ind
%             case 1 % Line 1-2
%                 p=linspace(P1,P2,n_pt);
%                 h=zeros(size(p));
%                 s=zeros(size(p));
%                 s(1)=S1;
%                 h(1)=H1;
%                 for k=2:n_pt-1
%                     [~,h(k),s(k),~,~] = compression(options.eta_SiC,p(k),H1,S1);
%                 end
%                 h(end)=H2;
%                 s(end)=S2;
%                 plot(s,h,'k','HandleVisibility','off');
%                 
%             case 2 % Line 2-8'
%                 p=linspace(P2,P8_prime,n_pt);
%                 s=linspace(S2,S8_prime,n_pt);
%                 h=zeros(size(s));
%                 h(1)=H2;
%                 for k=2:n_pt-1
%                     h(k)=XSteam('h_ps',p(k),s(k));
%                 end
%                 h(end)=H8_prime;
%                 plot(s,h,'k','HandleVisibility','off');
%                 
%             case 3 % Line 8''-8
%                 p=linspace(P8_seconde,P8,n_pt);
%                 s=linspace(S8_seconde,S8,n_pt);
%                 h=zeros(size(s));
%                 h(1)=H8_seconde;
%                 for k=2:n_pt-1
%                     h(k)=XSteam('h_ps',p(k),s(k));
%                 end
%                 h(end)=H8;
%                 plot(s,h,'k','HandleVisibility','off');
%                 
%             case 4 % Line 8-6
%                 p=linspace(P8,P6,n_pt);
%                 s=linspace(S8,S6,n_pt);
%                 h=zeros(size(s));
%                 h(1)=H8;
%                 for k=2:n_pt-1
%                     h(k)=XSteam('h_ps',p(k),s(k));
%                 end
%                 h(end)=H6;
%                 plot(s,h,'k','HandleVisibility','off');
%                 
%             case 5 % Line 10''-3
%                 p=linspace(P10_seconde,P3,n_pt);
%                 s=linspace(S10_seconde,S3,n_pt);
%                 h=zeros(size(s));
%                 h(1)=H10_seconde;
%                 for k=2:n_pt-1
%                     h(k)=XSteam('h_ps',p(k),s(k));
%                 end
%                 h(end)=H3;
%                 plot(s,h,'k','HandleVisibility','off');
%                 
%             case 6 % Line 8'-9'
%                 p=linspace(P8_prime,P9_prime,n_pt);
%                 h=zeros(size(p));
%                 s=zeros(size(p));
%                 s(1)=S8_prime;
%                 h(1)=H8_prime;
%                 for k=2:n_pt-1
%                     [~,h(k),s(k),~,~] = compression(options.eta_SiC,p(k),H8_prime,S8_prime);
%                 end
%                 h(end)=H9_prime;
%                 s(end)=S9_prime;
%                 plot(s,h,'k','HandleVisibility','off');
%                 
%             case 7 % Line 9'-10'
%                 p=linspace(P9_prime,P10_prime,n_pt);
%                 h=zeros(size(p));
%                 s=zeros(size(p));
%                 s(1)=S9_prime;
%                 h(1)=H9_prime;
%                 for k=2:n_pt-1
%                     [~,h(k),s(k),~,~] = compression(options.eta_SiC,p(k),H9_prime,S9_prime);
%                 end
%                 h(end)=H10_prime;
%                 s(end)=S10_prime;
%                 plot(s,h,'k','HandleVisibility','off');
%                 
%             case 8 % Line 8'-8''
%                 p=linspace(P8_prime,P8_seconde,n_pt);
%                 s=linspace(S8_prime,S8_seconde,n_pt);
%                 h=zeros(size(s));
%                 h(1)=H8_prime;
%                 for k=2:n_pt-1
%                     h(k)=XSteam('h_ps',p(k),s(k));
%                 end
%                 h(end)=H8_seconde;
%                 plot(s,h,'k','HandleVisibility','off');
%                 
%             case 9 % Line 9'-9''
%                 p=linspace(P9_prime,P9_seconde,n_pt);
%                 s=linspace(S9_prime,S9_seconde,n_pt);
%                 h=zeros(size(s));
%                 h(1)=H9_prime;
%                 for k=2:n_pt-1
%                     h(k)=XSteam('h_ps',p(k),s(k));
%                 end
%                 h(end)=H9_seconde;
%                 plot(s,h,'k','HandleVisibility','off');
%                 
%             case 10 % Line 10'-10''
%                 p=linspace(P10_prime,P10_seconde,n_pt);
%                 s=linspace(S10_prime,S10_seconde,n_pt);
%                 h=zeros(size(s));
%                 h(1)=H10_prime;
%                 for k=2:n_pt-1
%                     h(k)=XSteam('h_ps',p(k),s(k));
%                 end
%                 h(end)=H10_seconde;
%                 plot(s,h,'k','HandleVisibility','off');
%                 
%             case 11 %3-4
%                 p=linspace(P3,P4,n_pt);
%                 h=zeros(size(p));
%                 s=zeros(size(p));
%                 s(1)=S3;
%                 h(1)=H3;
%                 for k=2:n_pt-1
%                     [~,h(k),s(k),~,~] = detente(options.eta_SiT,p(k),H3,S3);
%                 end
%                 h(end)=H4;
%                 s(end)=S4;
%                 plot(s,h,'k','HandleVisibility','off');
%                 
%             case 12 %5-6
%                 p=linspace(P5,P6,n_pt);
%                 h=zeros(size(p));
%                 s=zeros(size(p));
%                 s(1)=S5;
%                 h(1)=H5;
%                 for k=2:n_pt-1
%                     [~,h(k),s(k),~,~] = detente(options.eta_SiT,p(k),H5,S5);
%                 end
%                 h(end)=H6;
%                 s(end)=S6;
%                 plot(s,h,'k','HandleVisibility','off');
%             case 13 %6-7
%                 p=linspace(P6,P7,n_pt);
%                 h=zeros(size(p));
%                 s=zeros(size(p));
%                 s(1)=S6;
%                 h(1)=H6;
%                 for k=2:n_pt-1
%                     [~,h(k),s(k),~,~] = detente(options.eta_SiT,p(k),H6,S6);
%                 end
%                 h(end)=H7;
%                 s(end)=S7;
%                 plot(s,h,'k','HandleVisibility','off');
%             case 14 %9''-9
%                 p=linspace(P9_seconde,P9,n_pt);
%                 s=linspace(S9_seconde,S9,n_pt);
%                 h=zeros(size(s));
%                 h(1)=H9_seconde;
%                 for k=2:n_pt-1
%                     h(k)=XSteam('h_ps',p(k),s(k));
%                 end
%                 h(end)=H9;
%                 plot(s,h,'k','HandleVisibility','off');
%             case 15 %9-4
%                 p=linspace(P9,P4,n_pt);
%                 s=linspace(S9,S4,n_pt);
%                 h=zeros(size(s));
%                 h(1)=H9;
%                 for k=2:n_pt-1
%                     h(k)=XSteam('h_ps',p(k),s(k));
%                 end
%                 h(end)=H4;
%                 plot(s,h,'k','HandleVisibility','off');
%             case 16 %2-7
%                 p=linspace(P2,P7,n_pt);
%                 s=linspace(S2,S7,n_pt);
%                 h=zeros(size(s));
%                 h(1)=H2;
%                 for k=2:n_pt-1
%                     h(k)=XSteam('h_ps',p(k),s(k));
%                 end
%                 h(end)=H7;
%                 plot(s,h,'k','HandleVisibility','off');
%             case 17 %4-5
%                 p=linspace(P4,P5,n_pt);
%                 s=linspace(S4,S5,n_pt);
%                 h=zeros(size(s));
%                 h(1)=H4;
%                 for k=2:n_pt-1
%                     h(k)=XSteam('h_ps',p(k),s(k));
%                 end
%                 h(end)=H5;
%                 plot(s,h,'k','HandleVisibility','off');
%         end
%     end


end