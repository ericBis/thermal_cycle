%%%CHANGEMENT%%%
%r - kcc - kmeca - P_e - T_0 - T_3 - eta_PiC - eta_PiT
%T0 - T_STmax - eta_mec - pdrum - pmid - x7 - eta_SiC - eta_SiT
%Parametres à faire varier sont: pmid puis pdrum tel que pdrum < pmid =>
%influence sur le rdt ST exerg (tqt on ait un titre acceptable et déterminer à partir de quel titre c'est pas acceptable)
%Puis pour cette plage faire varier les pressions pour voir comment varier
%le rdt STex et totex 
%
%Pour varier la température du condenseur et ses conséquence sur les rdt ST
%exerg et totex et pie chart
%
%Faire varier les parametres (r,kcc,P_e, etaPiC,etaPit) de la GT qui influencerai la chaudiere.
%Faire 2 puis plus de paramètres et voir les conséquences sur la chaudière
%
function GivePlotCCGT()
variation_PDRUM=(1:0.5:10);
longeurVecteur=length(variation_PDRUM);
valeur_RdtExergST=variation_PDRUM;%preallocation
valeur_RdtEnergTotCCGT=variation_PDRUM;%preallocation
valeur_RdtEnergST=variation_PDRUM;%preallocation
valeur_RdtExergTotCCGT=variation_PDRUM;%preallocation
for i=1:longeurVecteur
    P_eg = 225e3; %[kW]
    display=0;
    options=struct;
    options.pdrum=variation_PDRUM(i); %[bar]:Drum pressure
    [ETA,~,~] = CCGT(P_eg,options,display);
    valeur_RdtEnergST(i)=ETA(1);
    valeur_RdtExergST(i)=ETA(4);
    valeur_RdtEnergTotCCGT(i)=ETA(3);
    valeur_RdtExergTotCCGT(i)=ETA(6);
end
figure()
hold on
plot(variation_PDRUM,valeur_RdtEnergST,valeur_RdtExergST);
xlabel('pdrum');
ylabel('rendement');
legend('rdt energ.','rdt exerg.');
title('CCGT');
end