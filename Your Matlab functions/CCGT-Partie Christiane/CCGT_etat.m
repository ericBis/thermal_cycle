%Valeurs initiales de la CCGT
  P_EG = 225e3; %[kW]
  OPTIONS=struct;
  OPTIONS.T0 =15; %[°C]:Reference temperature
  OPTIONS.T_ext=25; %[°C]:External temperature
  OPTIONS.T_STmax=565; %[°C]:maximum temperature on ST cycle
  OPTIONS.eta_mec=0.956; %[-]:mecanic efficiency of shafts bearings
  OPTIONS.pdrum=4; %[bar]:Drum pressure
  OPTIONS.pmid=28; %[bar]:Intermediary pressure level
  OPTIONS.x7=0.95;
  OPTIONS.eta_SiC=0.9;
  OPTIONS.eta_SiT=0.9;
  OPTIONS.GT=struct;%[struct]
             %%%%DESCRIPION DE LA STRUCTURE GT%%%%%
             OPTIONS.GT.k_mec=0.015;
             OPTIONS.GT.T_0=15;
             OPTIONS.GT.T_ext=15;
             OPTIONS.GT.r=15;
             OPTIONS.GT.k_cc=0.05;  
             OPTIONS.GT.T_3=1250;   
             OPTIONS.GT.eta_PiC=OPTIONS.eta_SiC;
             OPTIONS.GT.eta_PiT=OPTIONS.eta_SiT;
             OPTIONS.GT.T4G=615;%A RETIRER ET SAVOIR COMMENT COLLECTER T4g et Y ACCEDER
             DISPLAY = 1;   %If 1,plot. If 0,do not plot.
  
  %%%%Calcul des etats de la CCGT%%%%
  %Variable etat 0:
  H0=XSteam('hL_T',OPTIONS.T0);
  T0=273.15+OPTIONS.T0;%[K]
  S0=XSteam('sL_T',OPTIONS.T0);
  
  %Pinch et DeltaEchangeur
  DeltaTa=OPTIONS.GT.T4G-OPTIONS.T_STmax;%[°C] difference de temperature entre T4g et T3
  DeltaPinch=15;%[°C] difference de temperature entre fumees et eau co_courante
  DeltaCondenseur=21;%[°C] difference de temperature entre riviere et vapeur
  
  %%%Etats du cycle
  T3=OPTIONS.T_STmax;%[°C] difference fumee en sortie de GT et eau
  P3=110;%[bar] valeur imposee
  H3=XSteam('h_pT',P3,T3);%[kJ/kg]
  S3=XSteam('s_pT',P3,T3);%[kJ/kgK]
  E3=(H3-H0)-T0*(S3-S0);%[kJ/kg]
  
  T5=T3;%[°C] IMPOSEE
  P5=28;%[bar] IMPOSEE
  H5=XSteam('h_pT',P5,T5);
  S5=XSteam('s_pT',P5,T5);
  E5=(H5-H0)-T0*(S5-S0);
  
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
  P10_prime=P3;%[bar]
  P10_seconde=P3;%[bar]
  T10_prime=XSteam('Tsat_p',P10_prime);
  T10_seconde=T10_prime;
  X10_prime=0;
  X10_seconde=1;
  H10_prime=XSteam('hL_p',P10_prime);
  H10_seconde=XSteam('hV_p',P10_seconde);
  S10_prime=XSteam('sL_p',P10_prime);
  S10_seconde=XSteam('sV_p',P10_seconde);
  E10_prime=(H10_prime-H0)-T0*(S10_prime-S0) ;
  E10_seconde=(H10_seconde-H0)-T0*(S10_seconde-S0) ;
  
  P6=OPTIONS.pdrum;%[bar] IMPOSEE
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
  
  S6s=S5;
  P6s=P6;
  H6s=XSteam('h_ps',P6s,S6s);
  H6=H5-OPTIONS.eta_SiT*(H5-H6s);
  S6=XSteam('s_ph',P6,H6);
  T6=XSteam('T_ps',P6,S6);
  E6=(H6-H0)-T0*(S6-S0);
  
  P4=P5;
  P4s=P4;
  S4s=S3;
  H4s=XSteam('h_ps',P4s,S4s);
  H4=H3-OPTIONS.eta_SiT*(H3-H4s);
  S4=XSteam('s_ph',P4,H4);
  T4=XSteam('T_ps',P4,S4);
  E4=(H4-H0)-T0*(S4-S0);
  
  P9=P5;
  T9=T10_prime;
  H9=XSteam('h_pT',P9,T9);
  S9=XSteam('s_pT',P9,T9);
  E9=(H9-H0)-T0*(S9-S0);
  
  X7=OPTIONS.x7;
  T7=DeltaCondenseur+T0;
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
  H2=H1+(1/OPTIONS.eta_SiC)*(H2s-H1);
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
  
  %%%%MASSFLOW de la CCGT%%%%
  massflowWater_HP=1;
  massflowWater_IP=1;
  massflowWater_LP=1;
  x0=[massflowWater_HP,massflowWater_IP,massflowWater_LP];
  funFlow=@CalculationFlow;
  mHPmMIPmLP=fsolve(@(x)funFlow(x),x0);

%   %%%%%TO DO%%%%%%%%%%
%   -Recuperer MFumee, MFComb et MFair de la GT;
%   -Calculer h5g et donc T5g;
%   -Calculer tout les rendements;

% ETA=ones(1,11);
% % ETA is a vector with :
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
%MASSFLOW=ones(1,5);




  
  %Recuperer MFumee,et definir les enthalpies comme variables globale
  function systemMassFlow = CalculationFlow(MassflowIterInitiale)
  systemMassFlow=ones(3,1);
  systemMassFlow(1,1)= -MFumee*(h_4g-h_HPg)+MassflowIterInitiale(1)*(H3-H10_prime+H5-H4)+MassflowIterInitiale(2)*(H5-H9);
  systemMassFlow(2,1)= -MFumee*(h_HPg-h_IPg)+MassflowIterInitiale(1)*(H10_prime-H9_prime)+MassflowIterInitiale(2)*(H9-H9_prime)+MassflowIterInitiale(3)*(H6-H8);
  systemMassFlow(3,1)= -MFumee*(h_IPg-h_LPg)+MassflowIterInitiale(1)*(H9_prime-H8_prime)+MassflowIterInitiale(2)*(H9_prime-H8_prime)+MassflowIterInitiale(3)*(H8-H8_prime);
  end
  
  
  
  
  