%%%% VARIATION DE 1000 à 1600 degree pour voir l'influence de la chaudiere
%%%% avec le rapport de compression
function GivePlotCCGT4()
P_eg=225e3;
options=struct;
options.pmid=20;
options.pdrum=4;
options.x7=0.88;
options.GT=struct;
Vecteur_R=(12:0.5:30);
VecteurTemp=(1000:200:1600);
VecteurEtaExergTot=VecteurTemp;
display=0;


for i=1:length(VecteurTemp)
    options.GT.T_3=VecteurTemp(i);
    for j=1:length(Vecteur_R)
    options.GT.r=Vecteur_R(j);
    [ETA,~,~] = CCGT(P_eg,options,display);
    VecteurEtaExergTot(j)=ETA(6);
    end
        hold on
        plot(Vecteur_R,VecteurEtaExergTot);
        title('CCGT: Rdt. totex en fonction de r')
        xlabel('r')
        ylabel('Rdt. totex')
end
VecteurEtaExergTot(end)%0.6007
legend('t3GT=1000°C','t3GT=1200°C','t3GT=1400°C','t3GT=1600°C');
end