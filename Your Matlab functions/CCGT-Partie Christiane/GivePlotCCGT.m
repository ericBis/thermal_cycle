%%%CHANGEMENT%%%
%Variation de pmid et pdrum tel que ensemble le rdt soit maximal
function GivePlotCCGT()
    P_eg = 225e3; %[kW]
    display=0;
    options=struct;
    variation_PDRUM=(1:4);
    variation_PMID=(20:0.5:28);
    valeur_RdtExergTotCCGT=variation_PMID;
    figure();
    for i=1:length(variation_PDRUM)
        options.pdrum=variation_PDRUM(i);
        for j=1:length(variation_PMID)
            options.pmid=variation_PMID(j);
            [ETA,~,~] = CCGT(P_eg,options,display);
            valeur_RdtExergTotCCGT(j)=ETA(6);
        end
        hold on
        plot(variation_PMID,valeur_RdtExergTotCCGT);
        xlabel('pmid [bar]');
        ylabel('Rdt. totex');
        title('CCGT: Rdt. totex en fonction de pmid');
    end
valeur_RdtExergTotCCGT(1) %rdt=0.5322 avec pmid=20 et pdrum=4
legend('pdrum=1[bar]','pdrum=2[bar]','pdrum=3[bar]','pdrum=4[bar]');
end