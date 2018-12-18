%%%CHANGEMENT
% Variation de R et de eta_PIC=eta-PIT pour voir le rendement exerg tot
%On voit que pour 0.92 et r=30 ou est bon

function GivePlotCCGT2()
P_eg=225e3;
options=struct;
options.pdrum=4;
options.pmid=20;
options.x7=0.88;
options.GT=struct;
Vecteur_R=(12:0.5:30);
Vecteur_Rdt_ExergTot=Vecteur_R;
Vecteur_Eta_Pi=(0.8:0.04:0.92);
figure()
    for j=1:length(Vecteur_Eta_Pi)
        options.GT.eta_PiC=Vecteur_Eta_Pi(j);
        options.GT.eta_PiT=Vecteur_Eta_Pi(j);
        for i=1:length(Vecteur_R)
        display=0;
        options.GT.r=Vecteur_R(i);
        [ETA,~,~] = CCGT(P_eg,options,display);
        Vecteur_Rdt_ExergTot(i)=ETA(6);
        end
        hold on
        plot(Vecteur_R,Vecteur_Rdt_ExergTot);
        title('CCGT: Rdt. totex en fonction de r')
        xlabel('r')
        ylabel('Rdt. totex')
    end
    Vecteur_Rdt_ExergTot(end) %0.5666
    legend('eta_PiC=eta_PiT=0.80','eta_PiC=eta_PiT=0.84','eta_PiC=eta_PiT=0.88','eta_PiC=eta_PiT=0.92');
end
