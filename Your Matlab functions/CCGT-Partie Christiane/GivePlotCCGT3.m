%%%% CHANGEMENT FAIRE VARIER le x7 de 88 à 95 pour pmid=20bar et pdrum=4
function GivePlotCCGT3()
P_eg=225e3;
options=struct;
options.pmid=20;
options.pdrum=4;
display=0;
n=20;
VecteurTitre=(0.88:(0.95-0.88)/n:0.95);%Vecteur à 5 valeurs
VecteurEtaExergTot=VecteurTitre;

for i=1:length(VecteurTitre)
options.x7=VecteurTitre(i);
[ETA,~,~] = CCGT(P_eg,options,display);
VecteurEtaExergTot(i)=ETA(6);
end
VecteurEtaExergTot(1)%0.5551
figure()
plot(VecteurTitre,VecteurEtaExergTot)
title('CCGT: Rdt. totex en fonction de x7')
xlabel('x7[-]')
ylabel('Rdt. totex')
end