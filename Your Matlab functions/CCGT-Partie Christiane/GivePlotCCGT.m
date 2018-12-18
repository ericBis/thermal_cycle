%%%CHANGEMENT%%%
%r - kcc - kmeca - P_e - T_0 - T_3 - eta_PiC - eta_PiT
%T0 - T_STmax - eta_mec - pdrum - pmid - x7 - eta_SiC - eta_SiT
%
%Pour varier la température du condenseur et ses conséquence sur les rdt ST
%exerg et totex et pie chart
%
%Faire varier les parametres (r,kcc,P_e, etaPiC,etaPit) de la GT qui influencerai la chaudiere.
%Faire 2 puis plus de paramètres et voir les conséquences sur la chaudière
%
function GivePlotCCGT()

%VARIATION DU DRUM
variation_PDRUM=(1:0.5:4);
longeurVecteur=length(variation_PDRUM);
valeur_RdtEnergTotCCGT=variation_PDRUM;%preallocation
valeur_RdtExergTotCCGT=variation_PDRUM;%preallocation

for i=1:longeurVecteur
    P_eg = 225e3; %[kW]
    display=0;
    options=struct;
    options.pdrum=variation_PDRUM(i); %[bar]:Drum pressure
    [ETA,~,~] = CCGT(P_eg,options,display);
    valeur_RdtEnergTotCCGT(i)=ETA(3);
    valeur_RdtExergTotCCGT(i)=ETA(6);
end
RDT_ENERG_MAX_PDRUM_CCGT=valeur_RdtEnergTotCCGT(end);%rdt quand pdrum=4
RDT_EXERG_MAX_PDRUM_CCGT=valeur_RdtExergTotCCGT(end);%rdt quand pdrum=4

figure(1)
subplot(1,2,1);
plot(variation_PDRUM,valeur_RdtEnergTotCCGT,'--',variation_PDRUM,valeur_RdtExergTotCCGT,':');
xlabel('pdrum');
ylabel('rendement');
legend('rdt energ. CCGT','rdt exerg.CCGT');
title('CCGT');


%VARIATION DU MID
variation_PMID=(20:0.5:28);
longeurVecteur=length(variation_PMID);
valeur_RdtEnergTotCCGT=variation_PMID;%preallocation
valeur_RdtExergTotCCGT=variation_PMID;%preallocation

for i=1:longeurVecteur
    P_eg = 225e3; %[kW]
    display=0;
    options=struct;
    options.pmid=variation_PMID(i); %[bar]:Drum pressure
    [ETA,~,~] = CCGT(P_eg,options,display);
    valeur_RdtEnergTotCCGT(i)=ETA(3);
    valeur_RdtExergTotCCGT(i)=ETA(6);
end
RDT_ENERG_MAX_PDRUM_CCGT;
RDT_EXERG_MAX_PDRUM_CCGT;
RDT_ENERG_MAX_PMID_CCGT=valeur_RdtEnergTotCCGT(1);
RDT_EXERG_MAX_PMID_CCGT=valeur_RdtExergTotCCGT(1);
RDT_ENERG_MIN_PMID_CCGT=valeur_RdtEnergTotCCGT(end);
RDT_EXERG_MIN_PMID_CCGT=valeur_RdtExergTotCCGT(end);
subplot(1,2,2);
plot(variation_PMID,valeur_RdtEnergTotCCGT,'--',variation_PMID,valeur_RdtExergTotCCGT,':');
xlabel('pmid');
ylabel('rendement');
legend('rdt energ. CCGT','rdt exerg. CCGT');
title('CCGT');
% %%%Pour une valeur de pdrum=4bar et pmid=28bar
% RDT_ENERG_MAX_PDRUM_CCGT=0.5445
% RDT_EXERG_MAX_PDRUM_CCGT=0.5230


% %%%Pour une valeur de pmid=20bar et pdrum=4
% RDT_ENERG_MAX_PDRUM_CCGT=0.5541
% RDT_EXERG_MAX_PDRUM_CCGT=0.5322
% si pmid=28bar
% RDT_ENERG_MIN_PMID_CCGT=0.5445
% RDT_EXERG_MIN_PMID_CCGT=0.5230

% %%%Pour une valeur de pdrum=4bar et pmid=20bar
%     valeur_RdtEnergTotCCGT=0.5541
%     valeur_RdtExergTotCCGT=0.5322 => augmenter pdrum et reduire pmid

%VARIATION POUR DES VALEUR PRDUM ET PMID QUI MAXIMISE RDT_TOT_EX
    P_eg = 225e3; %[kW]
    display=0;
    options=struct;
    options.pdrum=4;
    options.pmid=20; %[bar]:Drum pressure
    [ETA,~,~] = CCGT(P_eg,options,display);
    valeur_RdtEnergTotCCGT=ETA(3)
    valeur_RdtExergTotCCGT=ETA(6)
    

    
end