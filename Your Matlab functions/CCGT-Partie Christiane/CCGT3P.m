function [ETA MASSFLOW] = CCGT3P(P_eg,options,display)
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
%   -options.T0       [°C] : Reference temperature
%   -options.T_ext    [°C] : External temperature
%   -options.T_STmax  [°C] : maximum temperature on ST cycle
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
  P_eg = 225e3; %[kW]
  options=struct;
  options.T0 =15; %[°C]:Reference temperature
  options.T_ext=15; %[°C]:External temperature
  options.T_STmax=565; %[°C]:maximum temperature on ST cycle
  options.eta_mec=0.96; %[-]:mecanic efficiency of shafts bearings
  options.pdrum=4; %[bar]:Drum pressure
  options.pmid=28; %[bar]:Intermediary pressure level
  options.x7=0.95;
  options.eta_SiC=0.857;
  options.eta_SiT=0.857;
  options.GT=struct;%[struct]
           options.GT.k_mec=0.015;
           options.GT.T_0=15;
           options.GT.T_ext=15;
           options.GT.r=15;
           options.GT.k_cc=0.95;  
           options.GT.T_3=1250;   
           options.GT.eta_PiC=0.9;
           options.GT.eta_PiT=0.9;
   end
end

if isfield(options,'T_0')
    T_0 = options.T_0;
else
    T_0 = 15; % [°C]
end

%%%%OUTPUT GT%%%%
  [ETA, DATEN, DATEX, DAT, MASSFLOW, COMBUSTION] = GT(P_eg,options.GT,display);
  ETAGT=ETA;
  MASSFLOWGT=MASSFLOW;

%%%%OUTPUT VARIABLE CCGT%%%%
  ETA=ones(1,11);
  MASSFLOW=ones(1,5);
  MASSFLOW(4)=MASSFLOWGT(1);%MFAir
  MASSFLOW(5)=MASSFLOWGT(2);%MFComb
  MFumee=MASSFLOWGT(3);%MFumee
  
%%%%FUMEE EN 4G%%%%
  T_4g=DAT(1,4);%[°C]
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
  DeltaTa=T_4g-options.T_STmax;%[°C] difference de temperature entre T4g et T3
  DeltaPinch=8;%[°C] difference de temperature entre fumees et eau co_courante
  DeltaCondenseur=21;%[°C] difference de temperature entre riviere et vapeur
  
  %%%Etats du cycle
  T3=options.T_STmax;%[°C] difference fumee en sortie de GT et eau
  
  %L'etat 5 est entierement connu sur base des inputs
  T5=T3;%[°C] IMPOSEE
  P5=options.pmid;%[bar]
  H5=XSteam('h_pT',P5,T5);
  S5=XSteam('s_pT',P5,T5);
  E5=(H5-H0)-T0*(S5-S0);
  
  %L'etat 9 est determine à partir de l'etat 5
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
  T_5g=h_5g/CpMoyen_4g;%[°C]
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

  %%%%PIE CHART EXERG%%%%
  figure(1)
  dataExerg = [IrrevExergCombustion,IrrevExergTransfThermique,PerteExergCheminee,IrrevComplexeRotoriqueCCGT,PerteExergCondenseur,PerteMecaCCGT,PuissanceEffectiveST,PuissanceEffectiveGT]*(10^(-3));
  labelsExerg = {strcat('Irrev. comb.:',num2str(dataExerg(1)),'MW'),strcat('Irrev. transf.:',num2str(dataExerg(2)),'MW'),strcat('Perte cheminee:',num2str(dataExerg(3)),'MW'),strcat('Irrev. rotor CCGT:',num2str(dataExerg(4)),'MW'),strcat('Perte condenseur:',num2str(dataExerg(5)),'MW'),strcat('Perte meca CCGT:',num2str(dataExerg(6)),'MW'),strcat('Puiss. TV:',num2str(dataExerg(7)),'MW'),strcat('Puiss. TG:',num2str(dataExerg(8)),'MW')};
  pie(dataExerg,labelsExerg);
  title(strcat('Puiss. exergétique primaire:',num2str(PuissanceExergPrimaire*(10^(-3))),'MW'));

  %%%%PIE CHART ENERG%%%%
  figure(2)
  dataEnerg = [PerteEnergCheminee,PerteEnergCondenseur,PerteMecaCCGT,PuissanceEffectiveST,PuissanceEffectiveGT]*(10^(-3));
  labelsEnerg = {strcat('Perte cheminee:',num2str(dataEnerg(1)),'MW'),strcat('Perte condenseur:',num2str(dataEnerg(2)),'MW'),strcat('Perte meca CCGT:',num2str(dataEnerg(3)),'MW'),strcat('Puiss. TV:',num2str(dataEnerg(4)),'MW'),strcat('Puiss. TG:',num2str(dataEnerg(5)),'MW')};
  pie(dataEnerg,labelsEnerg);
  title(strcat('Puiss. energétique primaire:',num2str(PuissanceEnergPrimaire*(10^(-3))),'MW'));
  
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
  NumerateurEtaTransexSG=(MASSFLOW(1)*E3-sum(MASSFLOW(1:3))*E2);
  DenominateurEtaTransexSG=MFumee*(E_4g-E_5g);
  ETA(11)=NumerateurEtaTransexSG/DenominateurEtaTransexSG;%eta_transex, Heat exchanger overall exergy efficiency
  ETA(8)=ETA(10)*ETA(11)* ETA(9); %eta_gex, Steam generator exergy efficiency
 
  function systemMassFlow = CalculationFlow(MassflowIterInitiale,MFumee,h_4g,h_HPg,h_IPg,h_LPg,Enthalpies)
  systemMassFlow=ones(3,1);
  systemMassFlow(1,1)= -MFumee*(h_4g-h_HPg)+MassflowIterInitiale(1)*(Enthalpies(3)-Enthalpies(14)+Enthalpies(5)-Enthalpies(4))+MassflowIterInitiale(2)*(Enthalpies(5)-Enthalpies(11));
  systemMassFlow(2,1)= -MFumee*(h_HPg-h_IPg)+MassflowIterInitiale(1)*(Enthalpies(14)-Enthalpies(12))+MassflowIterInitiale(2)*(Enthalpies(11)-Enthalpies(12))+MassflowIterInitiale(3)*(Enthalpies(6)-Enthalpies(8));
  systemMassFlow(3,1)= -MFumee*(h_IPg-h_LPg)+MassflowIterInitiale(1)*(Enthalpies(12)-Enthalpies(9))+MassflowIterInitiale(2)*(Enthalpies(12)-Enthalpies(9))+MassflowIterInitiale(3)*(Enthalpies(8)-Enthalpies(9));
  end

end