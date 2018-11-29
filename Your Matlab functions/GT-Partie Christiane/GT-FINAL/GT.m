function [ETA, DATEN, DATEX, DAT, MASSFLOW, COMBUSTION] = GT(P_e,options,display)
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
%   -options.T_0   [°C] : Reference temperature
%   -options.T_ext [°C] : External temperature
%   -options.r     [-] : Comperssion ratio
%   -options.k_cc  [-] : Coefficient of pressure losses due to combustion
%                        chamber
%   -options.T_3   [°C] : Temperature after combustion (before turbine)
%   -option.eta_PiC[-] : Intern polytropic efficiency (Rendement
%                        polytropique interne) for compression
%   -option.eta_PiT[-] : Intern polytropic efficiency (Rendement
%                        polytropique interne) for expansion
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
% dat = {T_1       , T_2       , T_3       , T_4; [°C]
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
           P_e=100;%100MW
       end
      %%%%%DESCRIPION DE LA STRUCTURE%%%%%
%       options = struct;
%       options.k_mec=0.95;
%       options.T_0=15;
%       options.T_ext=25; 
%       options.r=10;   
%       options.k_cc=0.015;  
%       options.T_3=1050;   
%       options.eta_PiC=0.9; 
%       options.eta_PiT=0.9; 
   end
end


% Exemple of how to use (isfield' to check if an option has been given (or
% not)
if isfield(options,'T_0')
    T_0 = options.T_0;
else
    T_0 = 15;   
end
    %%%%%%APPEL FONCTION ETAT%%%%%
    T0=273.15;
    T1=273.15+options.T_0;
    T3=273.15+options.T_3;
    fun0 = @calculNaT2;
    x00 = [1.2,295];
    solutionNaetT2 = lsqnonlin(fun0,x00); %Determination de na et T2
    %%%
    fun00 = @calculNgT4;
    xf =[700, 1.1];%x0(1)=T4,x0(2)=lambda
    X=0;
    Y=4;
    solutionT4etLambda = fsolve(@(x)fun00(x,Y,X),xf); %Determination de T4 et lambda
    lambda=solutionT4etLambda(2);
    %%%
    TemperatureEtats=[T1,solutionNaetT2(2),T3,solutionT4etLambda(1)];
    %%%
    p0=1;
    pressions = [p0,p0*options.r,p0*options.r*options.k_mec,p0];
    %%%
    enth = enthalpy(TemperatureEtats);%Determination du vect enthalpie
    entrop = entropy(TemperatureEtats,pressions);
    exergie = exergy(enth,entrop);
    %%%
    matriceEtats=ones(4,4);
    matriceEtats(1,:)=TemperatureEtats;
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
    eta_meca=options.k_mec;
    xiter=0.2;%[kg/s]
    funny = @massFlowRateAir;
    PCI=PouvoirCalorifiqueInf(Y,X);


    %Masse Flow Rate
    MFair = fsolve(@(x)funny(x,Pe,eta_meca,lambda,pouvoirComburivor,enthalpieVector),xiter);%give MFair en %[kgAir/s]
    MFcomb = massFlowCombustible(pouvoirComburivor,lambda,MFair);%give MFcomb en [kgComb/s]
    MFfumee = massFlowFumee(MFair,MFcomb);%give MFfumee en [kgFumee/s]
    MFComposition = massFlowCompositionsFumee(Y,X,lambda,MFcomb);%give VECTEUR=MFComposition en [kgComposants/s]

    %Puissance energetique
    Pexh = puissanceExhaust(lambda,pouvoirComburivor,enthalpieVector,MFair,MFfumee); %en [MW]
    Pm = puissanceMotrice(enthalpieVector,MFair,MFfumee); %en [MW]
    Pfm = puissanceFrottementMecanique(Pe,Pm)%en [MW]

    %Puissance exergetique
    Protex = puissanceRotorExerg(exergieVector,MFair,MFfumee) %en [MW]
    exergieComb = exergieCombustible (Y,X); %en [kJ/kg]
    Pcomb = puissanceCombExerg(MFcomb,exergieComb) %en [MW] EXERGY 
    PexhExerg = puissanceExhaustExerg(lambda,pouvoirComburivor,exergieVector,MFair,MFfumee) %en [MW]
    Pm_exerg = puissanceMexerg(exergieVector,MFair,MFfumee);%en [MW]
    PfmecaExerg = puissanceFrottementMecaExerg(Pe,Pm_exerg); %en [MW]

    %Rendement energetique et exergetique
    eta_totEn = rdtTotalEnergetique(Pe,MFcomb,PCI);
    eta_cyclEn = rdtCyclEnergetique(lambda,pouvoirComburivor,enthalpieVector);
    eta_cyclEx = rdtCyclExergetique(Pm,exergieVector,MFair,MFfumee);
    eta_totEx = rdtTotalExergetique(Pe,Pcomb);
    eta_rotEx = rdtRotorExergetique(Pm,exergieVector,MFair,MFfumee);
    eta_combEx = rdtCombExergetique(exergieVector,MFair,MFfumee,Pcomb);


    %%%%%OUTPUT VARIABLES%%%%%
    ETA=[eta_cyclEn,eta_totEn,eta_cyclEx,eta_totEx,eta_rotEx,eta_combEx];
    DATEN=[Pfm,Pexh]*1000;%en [kW]
    DATEX=[PfmecaExerg,Protex,Pcomb,PexhExerg]*1000;%en [kW]
    DAT=matriceEtats;
    MASSFLOW=[MFair,MFcomb,MFfumee];
    COMBUSTION = struct();
    COMBUSTION.LHV=PouvoirCalorifiqueInf(Y,X);
    COMBUSTION.ec=exergieCombustible(Y,X);
    COMBUSTION.lambda=lambda;
    cpFumeeAt400K=cpMoyenFumee(400,TemperatureEtats(4),Y,X,lambda)/(TemperatureEtats(4)-400); %kJ/kgK
    COMBUSTION.Cp_g=cpFumeeAt400K;
    COMBUSTION.fum=MFComposition;
    COMBUSTION;
    
    %%%% PLOTS %%%%%
    
    %Plot T-s et h_s (A FAIRE)
    
    %Pie chart
    data = [Pfm,Protex,Pcomb,PexhExerg];
    labels = {strcat('Mechanical losses:',num2str(data(1)),'MW'),strcat('Mechanical exergetic power:',num2str(data(2)),'MW'),strcat('Primary exergy input:',num2str(data(3)),'MW'),strcat('Exhaust exergy losses:',num2str(data(4)),'MW')};
    pie(data,labels);
    
    %%%%%%%%%%%%%%%%%%%%%%%ALL FUNCTIONS OF GT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%PARTIE ENTHALPIE/ENTROPIE/EXERGIE%%%%%
    function enth = enthalpy(T)
    %TemperatureEtats(vecteur LIGNE avec T1,T2,T3,T4) EN KELVIN
    %cpMoyen est un vecteur de cpMoyen calculer correctement pour chaque etat
    %et calculer à chaque fois par rapport à la référence T0
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
    k_cc=options.k_mec;
    Tcomb=T0+options.T_ext;
    r=options.r;
    T2=solutionNaetT2(2);
    T3=273.15+options.T_3;
    eqq=ones(2,1);
    cpMoyenComb=(35.26)/(MasseMolaireCombustible(Y,X));
    eqq(1,1)=log(x(1)/T3) - (eta_PiT*constanteRfumee(Y,X,x(2))/cpMoyenFumee(x(1),T3,Y,X,x(2)))*log(1/(k_cc*r));
    eqq(2,1)=(T3-T0)-((1-options.k_cc)*PouvoirCalorifiqueInf(Y,X)+x(2)*pouvoirComburivore(Y,X)*cpMoyenAir(T0,T2)*(T2-T0)+cpMoyenComb*(Tcomb-T0))/(((pouvoirComburivore(Y,X))*x(2)+1)*cpMoyenFumee(T0,T3,Y,X,x(2)));
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
    function MFair = massFlowRateAir(x,Pe,eta_meca,lambda,pouvoirComburivor,enthalpieVector)
    MFair=((Pe*10^(3))/eta_meca) - x*(1+(1/(lambda*pouvoirComburivor)))*(enthalpieVector(3)-enthalpieVector(4))+x*(enthalpieVector(2)-enthalpieVector(1));
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
    function Pexh = puissanceExhaust(lambda,pouvoirComburivor,enthalpieVector,MFair,MFfumee)
    Pexh=(MFfumee*(1+1/(lambda*pouvoirComburivor))*enthalpieVector(4)-MFair*enthalpieVector(1))*10^(-3);
    end

    function Pm = puissanceMotrice(enthalpieVector,MFair,MFfumee)
    Pm=(MFfumee*(enthalpieVector(3)-enthalpieVector(4)) - MFair*(enthalpieVector(2)-enthalpieVector(1)))*10^(-3);
    end

    function Pfm = puissanceFrottementMecanique(Pe,Pm)
    Pfm=Pm-Pe;
    end

    %Functions Exergetic Power 
    function Protex = puissanceRotorExerg(exergieVector,MFair,MFfumee)
    Protex=(MFfumee*(exergieVector(3)-exergieVector(4))-MFair*(exergieVector(2)-exergieVector(1)))*10^(-3);
    end

    function exergieComb = exergieCombustible (Y,X)
    if (X==0)&&(Y==4)
        exergieComb=52215;
    end
    end

    function Pcomb = puissanceCombExerg(MFcomb,exergieCombustible)%en [MW]
    Pcomb = (MFcomb*exergieCombustible)*10^(-3);
    end

    function PexhExerg = puissanceExhaustExerg(lambda,pouvoirComburivor,exergieVector,MFair,MFfumee)
    PexhExerg =(MFfumee*(1+1/(lambda*pouvoirComburivor))*exergieVector(4)-MFair*exergieVector(1))*10^(-3);
    end %A DEMANDER (verification formule)

    function Pm_exerg = puissanceMexerg(exergieVector,MFair,MFfumee)
    Pm_exerg=(MFfumee*(exergieVector(3)-exergieVector(4))-MFair*(exergieVector(2)-exergieVector(1)))*10^(-3);
    end

    function PfmecaExerg = puissanceFrottementMecaExerg(Pe,Pm_exerg)
    PfmecaExerg=Pm_exerg-Pe;
    end

    %Exergetic and Energetic fficiencies 
    function eta_totEn = rdtTotalEnergetique(Pe,MFcomb,PCI)
    eta_totEn=Pe/((MFcomb*PCI)*10^(-3));
    end

    function eta_cyclEn = rdtCyclEnergetique(lambda,pouvoirComburivor,enthalpieVector)
    eta_cyclEn=((1+1/(lambda*pouvoirComburivor))*(enthalpieVector(3)-enthalpieVector(4))- (enthalpieVector(2)-enthalpieVector(1)))/((1+1/(lambda*pouvoirComburivor))*enthalpieVector(3)-enthalpieVector(2));
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

    function eta_combEx = rdtCombExergetique(exergieVector,MFair,MFfumee,Pcomb)
    eta_combEx=(MFfumee*exergieVector(3)-MFair*exergieVector(2))*(10^(-3))/(Pcomb);
    end
end
