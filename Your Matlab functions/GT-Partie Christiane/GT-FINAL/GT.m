function [ETA, DATEN, DATEX, DAT, MASSFLOW, COMBUSTION, FIG] = GT(P_e,options,display)
% GT Gas turbine modelisation
% GT(P_e,options,display) compute the thermodynamics states for a Gas
% turbine based on several inputs (given in OPTION) and based on a given
% electricity production P_e. It returns the main results. It can as well
% plots graphs if input argument DISPLAY = true (<=> DISPLAY=1)
%
% INPUTS (some inputs can be dependent on others => only one of these 2 can
%         be activated)
% P_E = electrical power output target [kW]
% OPTIONS is a structure containing :
%   -options.k_mec [-] : Shaft losses
%   -options.T_0   [�C] : Reference temperature
%   -options.T_ext [�C] : External temperature
%   -options.r     [-] : Comperssion ratio
%   -options.k_cc  [-] : Coefficient of pressure losses due to combustion
%                        chamber
%   -options.T_3   [�C] : Temperature after combustion (before turbine)
%   -options.eta_PiC[-] : Intern polytropic efficiency (Rendement
%                        polytropique interne) for compression
%   -options.eta_PiT[-] : Intern polytropic efficiency (Rendement
%                        polytropique interne) for expansion
%   -options.NTU
%DISPLAY = 1 or 0. If 1, then the code should plot graphics. If 0, then the
%          do not plot.
%
%OUPUTS :
% ETA is a vector with :
%   -eta(1) : eta_cyclen, cycle energy efficiency
%   -eta(2) : eta_toten, overall energy efficiency
%   -eta(3) : eta_cyclex, cycle exegy efficiency
%   -eta(4) : eta_totex, overall exergie efficiency
%   -eta(5) : eta_rotex, compressor-turbine exergy efficiency
%   -eta(6) : eta_combex, Combustion exergy efficiency
%   FYI : eta(i) \in [0;1] [-]
% DATEN is a vector with :
%   -daten(1) : perte_mec [kW]
%   -daten(2) : perte_ech [kW]
% DATEX is a vector with :
%   -datex(1) : perte_mec [kW]
%   -datex(2) : perte_rotex [kW]
%   -datex(3) : perte_combex [kW]
%   -datex(4) : perte_echex  [kW]
% DAT is a matrix containing :
% dat = {T_1       , T_2       , T_3       , T_4; [�C]
%        p_1       , p_2       , p_3       , p_4; [bar]
%        h_1       , h_2       , h_3       , h_4; [kJ/kg]
%        s_1       , s_2       , s_3       , s_4; [kJ/kg/K]
%        e_1       , e_2       , e_3       , e_4;};[kJ/kg]
% MASSFLOW is a vector containing :
%   -massflow(1) = m_a, air massflow [kg/s]
%   -massflow(2) = m_c, combustible massflow [kg/s]
%   -massflow(3) = m_f, exhaust gas massflow [kg/s]
%
% COMBUSTION is a structure with :
%   -combustion.LHV    : the Lower Heat Value of the fuel [kJ/kg]
%   -combustion.e_c    : the combuistible exergie         [kJ/kg]
%   -combustion.er     : exergie sensible                 [kJ/kg]
%   -combustion.lambda : the air excess                   [-]
%   -combustion.Cp_g   : heat capacity of exhaust gas at 400 K [kJ/kg/K]
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
%


%% Your Work

% Exemple of how to use 'nargin' to check your number of inputs
if nargin<3
    display=1;
    if nargin<2
        options = struct();
        if nargin<1
            P_e=230e3;%100MW
        end
        %%%%%DESCRIPION DE LA STRUCTURE%%%%%
        options = struct;
        options.k_mec=0.015;
        options.T_0=15;
        options.T_ext=15;
        options.r=18;
        options.k_cc=0.95;
        options.T_3=1400;
        options.eta_PiC=0.9;
        options.eta_PiT=0.9;
        options.NTU=0;
        options.ER=0;
    end
end


% Exemple of how to use (isfield' to check if an option has been given (or
% not)
if isfield(options,'T_0')
    T_0 = options.T_0;
else
    T_0 = 15;
end

if ~isfield(options,'k_mec')
    options.k_mec=0.015;
end

if ~isfield(options,'T_ext')
    options.T_ext=15;
end

if ~isfield(options,'r')
    options.r=18;
end

if ~isfield(options,'k_cc')
    options.k_c=0.95;
end

if ~isfield(options,'T_3')
    options.T_3=1400;
end

if ~isfield(options,'eta_PiC')
    options.eta_PiC=0.9;
end

if ~isfield(options,'eta_PiT')
    options.eta_PiT=0.9;
end

if ~isfield(options,'NTU')
    options.NTU=4;
end

if ~isfield(options,'GraphCompute')
    options.GraphCompute=1;
end


%%%%%%APPEL FONCTION ETAT%%%%%
T0=273.15;
T1=273.15+T_0;
T3=273.15+options.T_3;
fun0 = @calculNaT2;
x00 = [1.2,295];
solutionNaetT2 = lsqnonlin(fun0,x00); %Determination de na et T2
%%%
fun00 = @calculNgT4;
xf =[700, 1.1];%x0(1)=T4,x0(2)=lambda
X=0;
Y=4;
cpMoyenComb=(35.26)/(MasseMolaireCombustible(Y,X));
solutionT4etLambda = fsolve(@(x)fun00(x,Y,X),xf); %Determination de T4 et lambda
lambda=solutionT4etLambda(2);
%%%
TemperatureEtats=[T1,solutionNaetT2(2),T3,solutionT4etLambda(1)];%en Kelvin
%%%
p0=1;
pressions = [p0,p0*options.r,p0*options.r*options.k_cc,p0];%en bar
%%%
enth = enthalpy(TemperatureEtats);%Determination du vect enthalpie kJ/kg
entrop = entropy(TemperatureEtats,pressions);%kJ/kgK
exergie = exergy(enth,entrop);%kJ/kg
%%%
matriceEtats=ones(5,4);
matriceEtats(1,:)=TemperatureEtats-273;
matriceEtats(2,:)=pressions;
matriceEtats(3,:)=enth;
matriceEtats(4,:)=entrop;
matriceEtats(5,:)=exergie;
%%%

%%%%%APPEL FONCTION MASSFLOWRATE/PUISSANCE/RENDEMENT%%%%%
Pe=(P_e)*10^(-3);%en [MW]
enthalpieVector=enth;
exergieVector=exergie;
pouvoirComburivor= pouvoirComburivore(Y,X);
PCI=PouvoirCalorifiqueInf(Y,X);


%Masse Flow Rate
MFair = massFlowRateAir(Pe,lambda,pouvoirComburivor,enthalpieVector);%give MFair en %[kgAir/s]
MFcomb = massFlowCombustible(pouvoirComburivor,lambda,MFair);%give MFcomb en [kgComb/s]
MFfumee = massFlowFumee(MFair,MFcomb);%give MFfumee en [kgFumee/s]
MFComposition = massFlowCompositionsFumee(Y,X,lambda,MFcomb);%give VECTEUR=MFComposition en [kgComposants/s]

%Puissance energetique
Pexh = puissanceExhaust(enthalpieVector,MFair,MFfumee); %en [MW]
Pm = puissanceMotrice(enthalpieVector,MFair,MFfumee); %en [MW]
Pfm = puissanceFrottementMecanique(enthalpieVector,MFair,MFfumee);%en [MW]

%Puissance exergetique
ea=exergie(2);
ecr=cpMoyenComb*((options.T_ext+T0)-T1)-cpMoyenComb*T1*log((options.T_ext+T0)/T1);
    if options.ER==1
    er=ea*((lambda*pouvoirComburivor)/(lambda*pouvoirComburivor+1))+ecr*(1/(lambda*pouvoirComburivor+1));
    else er=0;
    end
Protex = puissanceRotorExerg(MFair,MFfumee,entrop);
exergieComb = exergieCombustible (Y,X); %en [kJ/kg]
Pcomb = puissanceCombExerg(MFcomb,exergieComb); %en [MW] EXERGY
PexhExerg = puissanceExhaustExerg(exergieVector,MFair,MFfumee); %en [MW]
PerteExComb = perteAlaComb (MFfumee,MFcomb,exergieComb,exergieVector);
PfmecaExerg = Pfm; %en [MW]

%Rendement energetique et exergetique

eta_totEn = rdtTotalEnergetique(Pe,MFcomb,PCI);
eta_cyclEn = rdtCyclEnergetique(MFair,Pm,enthalpieVector,MFfumee);
eta_cyclEx = rdtCyclExergetique(Pm,exergieVector,MFair,MFfumee);
eta_totEx = rdtTotalExergetique(Pe,Pcomb);
eta_rotEx = rdtRotorExergetique(Pm,exergieVector,MFair,MFfumee);
eta_combEx = rdtCombExergetique(exergieVector,MFfumee,Pcomb);
eta_meca = rdtMecanique(Pe,Pm);

%%%%%OUTPUT VARIABLES%%%%%
ETA=[eta_cyclEn,eta_totEn,eta_cyclEx,eta_totEx,eta_rotEx,eta_combEx];
DATEN=[Pfm,Pexh]*1000;%en [kW]
DATEX=[PfmecaExerg,Protex,PerteExComb,PexhExerg]*1000;%en [kW]
DAT=matriceEtats;
MASSFLOW=[MFair,MFcomb,MFfumee];
COMBUSTION = struct();
COMBUSTION.LHV=PouvoirCalorifiqueInf(Y,X);
COMBUSTION.er=er;
COMBUSTION.e_c=exergieCombustible(Y,X);
COMBUSTION.lambda=lambda;
cpFumeeAt400K=cpMoyenFumee(400,TemperatureEtats(4),Y,X,lambda)/(TemperatureEtats(4)-400); %kJ/kgK
COMBUSTION.Cp_g=cpFumeeAt400K;
COMBUSTION.fum=MFComposition;

%%%% CALCUL ECHANGEUR DE CHALEUR %%%%
%Calcul de T2R en fonction du NTU

if options.NTU~=0
T2R=(TemperatureEtats(4)*options.NTU+TemperatureEtats(2))/(1+options.NTU);
t2R=T2R-273.15;
%Calcul de T5 en fonction d'un bilan d'energie
fun000=@calculTemperatureT5;
x000=85;%[K]
T5 = fsolve(@(x)fun000(x),x000);
deltaPinchT4_T2R=TemperatureEtats(4)-T2R%[K] deltaPinchT4_T2R < deltaT5_T2
deltaT5_T2=T5-TemperatureEtats(2)%[K]
eta_EchangeurCyclenGT = eta_ExchangeCyclenGT()
S5=entropy_plotFum(T5,pressions(4));
H5=cpMoyenFumee(T0,T5,Y,X,lambda)*(T5-T0);
end

%%%% PLOTS %%%%%
if options.GraphCompute==1
    
    if display ==1
        disp='on';
    else
        disp='off';
    end
    
    if options.NTU~=0
        n_points=5;
    else
        n_points=4;
    end
    
    FIG = gobjects(4,1);
    n_pt=30;
    t=zeros(n_points*n_pt,1);
    s=zeros(n_points*n_pt,1);
    h=zeros(n_points*n_pt,1);
    kn=1;
    
    for k=1:n_points
        [t(kn:kn+n_pt-1),s(kn:kn+n_pt-1),h(kn:kn+n_pt-1)]=plot_TS_TH(k,n_pt);
        kn=kn+n_pt;
    end
    
    %T-S
    FIG(1)=figure('visible',disp);
    hold on
    title('Gas turbine T-s diagram');
    ylabel('Temperature [�C]');
    xlabel('Entropy [kJ/kg/K]');
    leg{1}='';
    for k=1:4
        plot(entrop(k),TemperatureEtats(k),'*')
        leg{k}=['Etat ' num2str(k)];
    end
    if options.NTU~=0
        plot(S5,T5,'*')
        leg{5}=['Etat ' num2str(5)];
    end
    legend(leg)
    plot(s,t,'k','HandleVisibility','off');
    
    %H-S
    FIG(2)=figure('visible',disp);
    hold on
    title('Gas turbine h-s diagram');
    ylabel('Enthalpie [kJ/kg]');
    xlabel('Entropy [kJ/kg/K]');
    for k=1:4
        plot(entrop(k),enth(k),'*')
    end  
    if options.NTU~=0
        plot(S5,H5,'*')
    end
    legend(leg)
    plot(s,h,'k','HandleVisibility','off');
    
    %Pie Charts
    
    FIG(3)=figure('visible',disp);
    pie([P_e, DATEN],{sprintf('%s \n %.1f [MW]','Effective power',P_e*1e-3),...
        sprintf('%s \n %.1f [MW]','Mechanical losses',DATEN(1)*1e-3),...
        sprintf('%s \n %.1f [MW]','Exhaust loss',DATEN(2)*1e-3)});
    title(sprintf('Gas Turbine Energy pie chart \n Primary power :  %.1f [MW]',...
        MFcomb*COMBUSTION.LHV*1e-3));
    
    FIG(4)=figure('visible',disp);
    pie([P_e, DATEX],{sprintf('%s \n %.1f [MW]','Effective power',P_e*1e-3),...
        sprintf('%s \n %.1f [MW]','Mechanical losses',DATEX(1)*1e-3),...
        sprintf('%s \n %.1f [MW]','Turbine & compressor irreversibilities',DATEX(2)*1e-3),...
        sprintf('%s \n %.1f [MW]','Combustion irreversibility',DATEX(3)*1e-3)...
        sprintf('%s \n %.1f [MW]','Exhaust loss',DATEX(4)*1e-3)});
    title(sprintf('Gas Turbine Exergy pie chart \n Primary exergy flux :  %.1f [MW]',...
        MFcomb*COMBUSTION.e_c*1e-3));
   
end

%%%%%%%%%%%%%%%%%%%%%%%ALL FUNCTIONS OF GT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%PARTIE ENTHALPIE/ENTROPIE/EXERGIE%%%%%
    function enth = enthalpy(T)
        %TemperatureEtats(vecteur LIGNE avec T1,T2,T3,T4) EN KELVIN
        %cpMoyen est un vecteur de cpMoyen calculer correctement pour chaque etat
        %et calculer � chaque fois par rapport � la r�f�rence T0
        %(cpMoyenAir1,cpMoyenAir2,cpMoyenFumee3,cpMoyenFumee4) LIGNE
        cpmoyenAirEtFumee=ones(1,4);
        cpmoyenAirEtFumee(1)=cpMoyenAir(T0,T(1));
        cpmoyenAirEtFumee(2)=cpMoyenAir(T0,T(2));
        cpmoyenAirEtFumee(3)=cpMoyenFumee(T0,T(3),Y,X,lambda);
        cpmoyenAirEtFumee(4)=cpMoyenFumee(T0,T(4),Y,X,lambda);
        enth=cpmoyenAirEtFumee.*(T-T0);%rendre vecteur de meme taille que "TemperatureEtats"
    end

    function entrop = entropy(Temp,pressions)
        entrop=ones(1,4);
        fun1=@(T)CpAir(T)./T;
        fun2=@(p)constanteRair./p;
        fun3=@(T)CpFumee(T,Y,X,lambda)./T;
        fun4=@(p)constanteRfumee(Y,X,lambda)./p;
        
        entrop(1)=integral(fun1,T0,Temp(1)) - integral(fun2,p0,pressions(1));
        entrop(2)=integral(fun1,T0,Temp(2)) - integral(fun2,p0,pressions(2));
        entrop(3)=integral(fun3,T0,Temp(3)) - integral(fun4,p0,pressions(3));
        entrop(4)=integral(fun3,T0,Temp(4)) - integral(fun4,p0,pressions(4));
    end

    function exergie = exergy(h,s)
        h0=h(1);
        s0=s(1);
        exergie=ones(1,4);
        exergie(1)=(h(1)-h0)-T1*(s(1)-s0);
        exergie(2)=(h(2)-h0)-T1*(s(2)-s0);
        exergie(3)=(h(3)-h0)-T1*(s(3)-s0);
        exergie(4)=(h(4)-h0)-T1*(s(4)-s0);
    end

%%%%%PARTIE COMBUSTION%%%%%
    function eqq = calculNgT4(x,Y,X)
        eta_PiT=options.eta_PiT;%x(1)=T4 x(2)=lambda
        k_cc=options.k_cc;
        Tcomb=T0+options.T_ext;
        r=options.r;
        T2=solutionNaetT2(2);
        T3=273.15+options.T_3;
        eqq=ones(2,1);
        eqq(1,1)=log(x(1)/T3) - (eta_PiT*constanteRfumee(Y,X,x(2))/cpMoyenFumee(x(1),T3,Y,X,x(2)))*log(1/(k_cc*r));
        eqq(2,1)=(T3-T0)-((1-0.01)*PouvoirCalorifiqueInf(Y,X)+x(2)*pouvoirComburivore(Y,X)*cpMoyenAir(T0,T2)*(T2-T0)+cpMoyenComb*(Tcomb-T0))/(((pouvoirComburivore(Y,X))*x(2)+1)*cpMoyenFumee(T0,T3,Y,X,x(2)));
    end

    function cpMFumee = cpMoyenFumee(Ta,Tb,Y,X,lambda)%Valable pour les etats 2,3,4 et tient compte de la compostition de l'air
        fun=@CpFumee;
        cpMFumee=(integral(@(T)fun(T,Y,X,lambda),Ta,Tb))/(Tb-Ta);
    end

    function cpF = CpFumee(a,Y,X,lambda)%Valable pour les etats 3 et 4, utiliser janaf directement. Compositions de l'air:O2 N2, CO2, H2O
        T=a;%Longueur du vecteur a chaque appel
        vecteur=find(T<300);
        L=length(vecteur);
        vvecteur=find(T>=300);
        m=length(vvecteur);
        cpF=ones(1,length(T));
        
        MmO2=32e-3;%[kgO2/molO2]
        MmCO2=44e-3;
        MmN2=28e-3;
        MmH2O=18e-3;
        
        propO2=(lambda*(1+Y/4-X/2)-1+X/2-Y/4)/(4.76*lambda*(1+Y/4-X/2)+X/2+Y/4);%Propotion molO2/molFumee
        propCO2=1/(4.76*lambda*(1+Y/4-X/2)+X/2+Y/4);%Propotion molCO2/molFumee
        propH2O=(Y/2)/(4.76*lambda*(1+Y/4-X/2)+X/2+Y/4);%Propotion molH2O/molFumee
        propN2=(3.76*(lambda*(1+Y/4-X/2)))/(4.76*lambda*(1+Y/4-X/2)+X/2+Y/4);%Propotion molN2/molFumee
        
        concentrationMassiqueO2=propO2*MmO2/(propO2*MmO2+propN2*MmN2+propCO2*MmCO2+propH2O*MmH2O);%[kgO2/kgFumee]
        concentrationMassiqueN2=propN2*MmN2/(propO2*MmO2+propN2*MmN2+propCO2*MmCO2+propH2O*MmH2O);
        concentrationMassiqueCO2=propCO2*MmCO2/(propO2*MmO2+propN2*MmN2+propCO2*MmCO2+propH2O*MmH2O);
        concentrationMassiqueH2O=propH2O*MmH2O/(propO2*MmO2+propN2*MmN2+propCO2*MmCO2+propH2O*MmH2O);
        
        if L>0%on trouve des temperatures<300
            cpF(vecteur)=concentrationMassiqueO2*janaf('c','O2',300)+concentrationMassiqueN2*janaf('c','N2',300)+concentrationMassiqueCO2*janaf('c','CO2',300)+concentrationMassiqueH2O*janaf('c','H2O',300);%[kJ/kgFumeeK] %rempli le vecteur avec les bons indices quand T<300K
        end
        if m>0%executer quand T>=300K
            cpF(vvecteur)=concentrationMassiqueO2*janaf('c','O2',T(vvecteur))+concentrationMassiqueN2*janaf('c','N2',T(vvecteur))+concentrationMassiqueCO2*janaf('c','CO2',T(vvecteur))+concentrationMassiqueH2O*janaf('c','H2O',T(vvecteur));%[kJ/kgFumeeK] +concentrationMassiqueH2O*janaf('c','H2O',T(vvecteur))
        end
    end

    function RstarFumee = constanteRfumee(Y,X,lambda)
        R=8.314e-3;%[kJ/mol.K]
        MmO2=32e-3;%[kgO2/molO2]
        MmCO2=44e-3;
        MmN2=28e-3;
        MmH2O=18e-3;
        propO2=(lambda*(1+Y/4-X/2)-1+X/2-Y/4)/(4.76*lambda*(1+Y/4-X/2)+X/2-Y/4);%Propotion molO2/molFumee
        propCO2=1/(4.76*lambda*(1+Y/4-X/2)+X/2-Y/4);%Propotion molCO2/molFumee
        propH2O=(Y/2)/(4.76*lambda*(1+Y/4-X/2)+X/2+Y/4);%Propotion molH2O/molFumee
        propN2=(3.76*(lambda*(1+Y/4-X/2)))/(4.76*lambda*(1+Y/4-X/2)+X/2-Y/4);%Propotion molN2/molFumee
        RstarFumee=R/(propO2*MmO2+propN2*MmN2+propCO2*MmCO2+propH2O*MmH2O);%[kJ/kgFumeeK]
    end

    function pci = PouvoirCalorifiqueInf(Y,X)%tiree d'une formule empirique
        formuleEmpirique=393400+102250*Y-(X/(1+0.5*Y))*(111000+102250*Y);
        pci=formuleEmpirique/(MasseMolaireCombustible(Y,X));%donne le PCI en kJ/kgComb de type CHyOx c'est un scalaire
    end

    function ma1 = pouvoirComburivore(Y,X)%Sachant que X et Y sont les proprietes du fuel
        ma1=(1+(Y-2*X)/4)*(32+3.76*28.15)/(12+Y+16*X);%[kgair/kgcomb] c'est un scalaire
    end

    function MMComb = MasseMolaireCombustible(Y,X)%Pour un combustible de type CHyOx
        MmO=16;%[kg/kmolO]
        MmH=1;%[kg/kmolH]
        MmC=12;%[kg/kmolC]
        MMComb=MmC+Y*MmH+X*MmO;%[kgComb/kmolComb] c'est un scalaire
    end

%%%%%PARTIE ADMISSION AIR%%%%%
    function eq = calculNaT2(x)%eq est un vecteur colonne qui determine na et T2
        eta_PiC=options.eta_PiC;%x(1)=na et x(2)=T2
        r=options.r;
        eq=ones(2,length(x));
        eq(1,1)= (x(1)-1)/x(1) - (1/eta_PiC)*(constanteRair/cpMoyenAir(T1,x(2)));
        eq(2,1)= x(2)/T1 -r^((x(1)-1)/x(1));
    end

    function RstarAir = constanteRair()%RstarAir est un scalaire
        R=8.314e-3;%[kJ/mol.K]
        MmO2=32e-3;%[kg/mol]
        MmN2=28e-3;%[kg/mol]
        propO2=0.21;%Propotion en mole
        propN2=0.79;%Propotion en mole
        RstarAir=R/(MmO2*propO2+MmN2*propN2);%[kJ/kgAirK]
    end

    function cpMAir = cpMoyenAir(Ta,Tb)%cpMAir est un scalaire. Valable pour les etats 1 et 2 et tient compte de la compostition de l'air
        fun=@CpAir;
        cpMAir=(integral(fun,Ta,Tb))/(Tb-Ta);%[kJ/kgAirK]
    end

    function cpA = CpAir(a) %cpA est un vecteur de cp du melange pour tout temperature
        %Valable pour les etats 1 et 2. Compositions de l'air: O2 et N2, rend un VECTEUR
        T=a;%Longueur du vecteur a chaque appel
        x=find(T<300);
        L=length(x);
        xx=find(T>=300);
        m=length(xx);
        cpA=ones(1,length(T));
        if L>0%on trouve des temperatures<300
            cpA(x)=1.009;%rempli le vecteur avec les bons indices quand T<300K
        end
        if m>0%executer quand T>=300K
            MmO2=32e-3;%[kg/mol]
            MmN2=28e-3;%[kg/mol]
            propO2=0.21;%Propotion en mole
            propN2=0.79;%Propotion en mole
            concentrationMassiqueO2=(MmO2*propO2)/(MmN2*propN2+MmO2*propO2);%[kgO2/kgAir]
            concentrationMassiqueN2=(MmN2*propN2)/(MmN2*propN2+MmO2*propO2);%[kgN2/kgAir]
            cpA(xx)=concentrationMassiqueO2*janaf('c','O2',T(xx))+concentrationMassiqueN2*janaf('c','N2',T(xx));%[kJ/kgAirK]
        end
    end

%%%%%PARTIE MASSFLOWRATE/PUISSANCE/RENDEMENT%%%%%

%Functions Masse Flow Rate
    function MFair = massFlowRateAir(Pe,lambda,pouvoirComburivor,enthalpieVector)
        coef=lambda*pouvoirComburivor;
        MFair=(Pe*1e3)/((1+1/coef)*(enthalpieVector(3)-enthalpieVector(4))*(1-options.k_mec)-(enthalpieVector(2)-enthalpieVector(1))*(1+options.k_mec));
    end

    function MFcomb = massFlowCombustible(pouvoirComburivor,lambda,MFair)
        MFcomb=MFair/(lambda*pouvoirComburivor);
    end

    function MFfumee = massFlowFumee(MFair,MFcomb)
        MFfumee=MFair+MFcomb;
    end

    function MFComposition = massFlowCompositionsFumee(Y,X,lambda,MFcomb)
        w=lambda*(1+Y/4-X/2);%rapport du nombre de mole du compose sur le nombre de mole de CHYOX de la combustion
        coefO2=(w+X/2-Y/4-1);
        coefCO2=1;
        coefN2=3.76*w;
        coefH2O=Y/2;
        
        MmO2=32e-3;%[kgO2/molO2]
        MmCO2=44e-3;
        MmN2=28e-3;
        MmH2O=18e-3;
        
        MasseMolaireCombustible=16e-3;
        
        moleFlowCHYOX=MFcomb/MasseMolaireCombustible;
        moleFlowO2=moleFlowCHYOX*coefO2;
        moleFlowCO2=moleFlowCHYOX*coefCO2;
        moleFlowN2=moleFlowCHYOX*coefN2;
        moleFlowH2O=moleFlowCHYOX*coefH2O;
        
        MFComposition=ones(1,4);
        MFComposition(1)=MmO2*moleFlowO2;
        MFComposition(2)=MmN2*moleFlowN2;
        MFComposition(3)=MmCO2*moleFlowCO2;
        MFComposition(4)=MmH2O*moleFlowH2O;
    end

%Functions Energetic Power
    function Pexh = puissanceExhaust(enthalpieVector,MFair,MFfumee)
        Pexh=(MFfumee*enthalpieVector(4)-MFair*enthalpieVector(1))*10^(-3);
    end

    function Pm = puissanceMotrice(enthalpieVector,MFair,MFfumee)
        Pm=(MFfumee*(enthalpieVector(3)-enthalpieVector(4)) - MFair*(enthalpieVector(2)-enthalpieVector(1)))*10^(-3);
    end

    function Pfm = puissanceFrottementMecanique(enthalpieVector,MFair,MFfumee)
        Pfm=options.k_mec*((MFfumee*(enthalpieVector(3)-enthalpieVector(4))+MFair*(enthalpieVector(2)-enthalpieVector(1)))*10^(-3));
    end

%Functions Exergetic Power
    function Protex = puissanceRotorExerg(MFair,MFfumee,entrop)
    Protex=((MFfumee*T0*(entrop(4)-entrop(3)))+MFair*T0*(entrop(2)-entrop(1)))*10^(-3);
    end

    function exergieComb = exergieCombustible (Y,X)
        if (X==0)&&(Y==4)
            exergieComb=52215;
        end
    end

    function Pcomb = puissanceCombExerg(MFcomb,exergieCombustible)%en [MW]
        Pcomb = (MFcomb*exergieCombustible)*10^(-3);
    end

    function PexhExerg = puissanceExhaustExerg(exergieVector,MFair,MFfumee)
        PexhExerg =(MFfumee*exergieVector(4)-MFair*exergieVector(1))*10^(-3);
    end

 function PerteExComb = perteAlaComb (MFfumee,MFcomb,exergieComb,exergieVector)
    PerteExComb=(MFcomb*exergieComb-MFfumee*(exergieVector(3))+MFair*exergieVector(2))*10^(-3);
    end
%Exergetic and Energetic fficiencies
    function eta_totEn = rdtTotalEnergetique(Pe,MFcomb,PCI)
        eta_totEn=Pe/((MFcomb*PCI)*10^(-3));
    end

    function eta_cyclEn = rdtCyclEnergetique(MFair,Pm,enthalpieVector,MFfumee)
        Pc=(MFfumee*enthalpieVector(3)-MFair*enthalpieVector(2))*10^(-3);
        eta_cyclEn=Pm/Pc;
    end

    function eta_cyclEx = rdtCyclExergetique(Pm,exergieVector,MFair,MFfumee)
        eta_cyclEx=Pm/((MFfumee*(exergieVector(3))-MFair*(exergieVector(2)))*10^(-3));
    end

    function eta_totEx = rdtTotalExergetique(Pe,Pcomb)
        eta_totEx=Pe/Pcomb;
    end

    function eta_rotEx = rdtRotorExergetique(Pm,exergieVector,MFair,MFfumee)
        eta_rotEx=Pm/((MFfumee*(exergieVector(3)-exergieVector(4))-MFair*(exergieVector(2)-exergieVector(1)))*10^(-3));
    end

    function eta_combEx = rdtCombExergetique(exergieVector,MFfumee,Pcomb)
        eta_combEx=MFfumee*(exergieVector(3)-er)*(10^(-3))/(Pcomb);
    end

    function eta_meca = rdtMecanique(Pe,Pm)
        eta_meca=Pe/Pm;
    end

%%%%PARTIE Echangeur de chaleur%%%%

    function calculT5 = calculTemperatureT5(x)
        calculT5=MFair*cpMoyenAir(TemperatureEtats(2),T2R)*(T2R-TemperatureEtats(2))-MFfumee*cpMoyenFumee(TemperatureEtats(4),x,Y,X,lambda)*(TemperatureEtats(4)-x);
    end

    function eta_EchangeurCyclenGT = eta_ExchangeCyclenGT()
        coef1=(1+1/(lambda*pouvoirComburivor));
        coef2=cpMoyenFumee(TemperatureEtats(4),TemperatureEtats(3),Y,X,lambda)/cpMoyenAir(TemperatureEtats(1),TemperatureEtats(2));
        coef3=(TemperatureEtats(3)/TemperatureEtats(1))*(1-(TemperatureEtats(4)/TemperatureEtats(3)));
        coef4=(TemperatureEtats(2)/TemperatureEtats(1))-1;
        coef5=cpMoyenFumee(TemperatureEtats(2),TemperatureEtats(3),Y,X,lambda)/cpMoyenAir(TemperatureEtats(1),TemperatureEtats(2));
        coef6=TemperatureEtats(3)/TemperatureEtats(1);
        coef7=TemperatureEtats(2)/TemperatureEtats(1);
        coef8=cpMoyenAir(TemperatureEtats(2),T2R)/cpMoyenAir(TemperatureEtats(1),TemperatureEtats(2));
        coef9=options.NTU/(1+options.NTU);
        coef10=(TemperatureEtats(4)-TemperatureEtats(2))/TemperatureEtats(1);
        numerateur=coef1*coef2*coef3-coef4;
        denominateur=coef1*coef5*coef6-coef7-coef8*coef9*coef10;
        eta_EchangeurCyclenGT=numerateur/denominateur;
    end

    function [t,S,H]=plot_TS_TH(ind,n_pt)
        optis = optimset('Display','off');
        switch ind
            case 1 % Line 1-2
                p=linspace(pressions(1),pressions(2),n_pt);
                t=linspace(TemperatureEtats(1),TemperatureEtats(2),n_pt);
                S=zeros(size(p));
                H=zeros(size(p));
                S(1)=entrop(1);
                t(1)=TemperatureEtats(1);
                H(1)=enth(1);
                x0 = [1.2,295];
                for c=2:n_pt-1
                    NaT = lsqnonlin(@(x) calculNaT12_plot(p(c),x),...
                        x0,[0 0],[100 2000],optis); %Determination de na et T
                    t(c)=NaT(2);
                    x0 = NaT;
                    S(c)=entropy_plotAir(t(c),p(c));
                    H(c)=cpMoyenAir(T0,t(c))*(t(c)-T0);
                end
                t(end)=TemperatureEtats(2);
                S(end)=entrop(2);
                H(end)=enth(2);
                % plot(S,t,'k','HandleVisibility','off');
            case 2 % Line 2-3
                p=linspace(pressions(2),pressions(3),n_pt); % On suposse une chute de pression lineaire l'augmentation de temperature
                t=linspace(TemperatureEtats(2),TemperatureEtats(3),n_pt);
                S=zeros(size(p));
                H=zeros(size(p));
                S(1)=entrop(2);
                t(1)=TemperatureEtats(2);
                H(1)=enth(2);
                for c=2:n_pt-1
                    S(c)=entropy_plotFum(t(c),p(c));
                    H(c)=cpMoyenFumee(T0,t(c),Y,X,lambda)*(t(c)-T0);
                end
                t(end)=TemperatureEtats(3);
                S(end)=entrop(3);
                H(end)=enth(3);
                % plot(S,t,'k','HandleVisibility','off');
            case 3 % Line 3-4
                p=linspace(pressions(3),pressions(4),n_pt);
                t=linspace(TemperatureEtats(3),TemperatureEtats(4),n_pt);
                S=zeros(size(p));
                H=zeros(size(p));
                S(1)=entrop(3);
                t(1)=TemperatureEtats(3);
                H(1)=enth(3);
                x0 = [700, 1.1];
                for c=2:n_pt-1
                    NgT = lsqnonlin(@(x) calculNgT34_plot(p(c),x),x0,...
                        [0 0],[2000 5],optis); %Determination de na et T
                    t(c)=NgT(1);
                    x0 = NgT;
                    S(c)=entropy_plotFum(t(c),p(c));
                    H(c)=cpMoyenFumee(T0,t(c),Y,X,lambda)*(t(c)-T0);
                end
                t(end)=TemperatureEtats(4);
                S(end)=entrop(4);
                H(end)=enth(4);
                %   plot(S,t,'k','HandleVisibility','off');
            case 4 % Line 4-1 or 4-5
                if options.NTU~=0 % Heat Exchanger
                    p=linspace(pressions(4),pressions(4),n_pt); % On suppose l'�change de chaleur isobare 
                    t=linspace(TemperatureEtats(4),T5,n_pt);
                    S=zeros(size(p));
                    H=zeros(size(p));
                    S(1)=entrop(4);
                    t(1)=TemperatureEtats(4);
                    H(1)=enth(4);
                    for c=2:n_pt
                        S(c)=entropy_plotFum(t(c),p(c));
                        H(c)=cpMoyenFumee(T0,t(c),Y,X,lambda)*(t(c)-T0);
                    end
                else % No heat exchanger
                    p=linspace(pressions(4),pressions(1),n_pt); % On suposse une chute de pression lineaire l'augmentation de temperature
                    t=linspace(TemperatureEtats(4),TemperatureEtats(1),n_pt);
                    S=zeros(size(p));
                    H=zeros(size(p));
                    S(1)=entrop(4);
                    t(1)=TemperatureEtats(4);
                    H(1)=enth(4);
                    for c=2:n_pt-1
                        S(c)=entropy_plotFum(t(c),p(c));
                        H(c)=cpMoyenFumee(T0,t(c),Y,X,lambda)*(t(c)-T0);
                    end
                    t(end)=TemperatureEtats(1);
                    S(end)=entrop(1);
                    H(end)=enth(1);
                end
            case 5 % Line 5-1
                if options.NTU~=0
                    p=linspace(pressions(4),pressions(1),n_pt); % On suppose l'�change de chaleur isobare
                    t=linspace(T5,TemperatureEtats(1),n_pt);
                    S=zeros(size(p));
                    H=zeros(size(p));
                    S(1)=S5;
                    t(1)=T5;
                    H(1)=H5;
                    for c=2:n_pt
                        S(c)=entropy_plotFum(t(c),p(c));
                        H(c)=cpMoyenFumee(T0,t(c),Y,X,lambda)*(t(c)-T0);
                    end                    
                end
        end
    end

    function eq = calculNaT12_plot(p,x)
        %eq est un vecteur colonne qui determine na et T
        %x(1)=na et x(2)=T
        r=p/pressions(1);
        eq=ones(2,1);
        eq(1,1)= (x(1)-1)/x(1) - (1/options.eta_PiC)*(constanteRair/cpMoyenAir(T1,x(2)));
        eq(2,1)= x(2)/T1 -r^((x(1)-1)/x(1));
    end

    function eq = calculNgT34_plot(p,x)
        %eq est un vecteur colonne qui determine ng et T
        %x(1)=T x(2)=lambda
        r=options.r;
        T2=solutionNaetT2(2);
        T3=273.15+options.T_3;
        Tcomb=T0+options.T_ext;
        eq=ones(2,1);
        eq(1)=log(x(1)/T3) - (options.eta_PiT*constanteRfumee(Y,X,x(2))/...
            cpMoyenFumee(x(1),T3,Y,X,x(2)))*log(p/pressions(3));
        eq(2)=(T3-T0)-((1-0.01)*PouvoirCalorifiqueInf(Y,X)+x(2)*pouvoirComburivore(Y,X)*cpMoyenAir(T0,T2)*(T2-T0)+cpMoyenComb*(Tcomb-T0))/(((pouvoirComburivore(Y,X))*x(2)+1)*cpMoyenFumee(T0,T3,Y,X,x(2)));
    end

    function entrop = entropy_plotAir(Temp,p)
        fun1=@(T)CpAir(T)./T;
        fun2=@(p)constanteRair./p;
        
        entrop=integral(fun1,T0,Temp) - integral(fun2,p0,p);
        
    end

    function entrop = entropy_plotFum(Temp,p)
        fun3=@(T)CpFumee(T,Y,X,lambda)./T;
        fun4=@(p)constanteRfumee(Y,X,lambda)./p;
        
        entrop=integral(fun3,T0,Temp) - integral(fun4,p0,p);
    end

end
