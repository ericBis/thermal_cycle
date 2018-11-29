%AMELIORATION A FAIRE:
%-ajouter un échangeur de chaleur : augmenter la temperature en sortie de combustion. Checker les
%rendements, et les puissances qu'on en retire.
%-faire varier pour l'instant le ratio 'r' pour voir son impact
%-faire varier la temperature du combustible avant la combustion
%-augmenter le nombre d'etages de compression (attention h2 est proche de
%h1) mais tenir compte aussi de l'augmentation de la pression qui
%"pourrait" impacter la temperature.

%Variables
Pe=230;%en [MW]
lambda=2.3;
enthalpieVector=[15.1,443.4,1681.8,731.5];
exergieVector=[0,403.2,1317.9,343.9];
pouvoirComburivore=16.999;
k_meca=0.015;
x0=0.2;%[kg/s]
fun = @massFlowRateAir;
PCI=50150;
X=0;
Y=4;

%Masse Flow Rate
MFair = fsolve(@(x)fun(x,Pe,k_meca,lambda,pouvoirComburivore,enthalpieVector),x0)%give MFair en %[kgAir/s]
MFcomb = massFlowCombustible(pouvoirComburivore,lambda,MFair)%give MFcomb en [kgComb/s]
MFfumee = massFlowFumee(MFair,MFcomb)%give MFfumee en [kgFumee/s]
MFComposition = massFlowCompositionsFumee(Y,X,lambda,MFcomb)%give VECTEUR=MFComposition en [kgComposants/s]

%Puissance energetique
Pexh = puissanceExhaust(lambda,pouvoirComburivore,enthalpieVector,MFair,MFfumee) %en [MW]
Pm = puissanceMotrice(enthalpieVector,MFair,MFfumee); %en [MW]
Pfm = puissanceFrottementMecanique(k_meca,enthalpieVector,MFair,MFfumee)%en [MW]

%Puissance exergetique
Protex = puissanceRotorExerg(exergieVector,MFair,MFfumee); %en [MW]
exergieComb = exergieCombustible (Y,X); %en [kJ/kg]
Pcomb = puissanceCombExerg(MFcomb,exergieComb) %en [MW] EXERGY 
PexhExerg = puissanceExhaustExerg(lambda,pouvoirComburivore,exergieVector,MFair,MFfumee) %en [MW]
Pm_exerg = puissanceMexerg(exergieVector,MFair,MFfumee);%en [MW]
PfmecaExerg = puissanceFrottementMecaExerg(Pe,Pm_exerg); %en [MW]

%Rendement energetique et exergetique
eta_totEn = rdtTotalEnergetique(Pe,MFcomb,PCI);
eta_cyclEn = rdtCyclEnergetique(lambda,pouvoirComburivore,enthalpieVector);
eta_cyclEx = rdtCyclExergetique(Pm,exergieVector,MFair,MFfumee);
eta_totEx = rdtTotalExergetique(Pe,Pcomb);
eta_rotEx = rdtRotorExergetique(Pm,exergieVector,MFair,MFfumee);
eta_combEx = rdtCombExergetique(exergieVector,MFair,MFfumee,Pcomb);
eta_meca = rdtMeca (Pfm,Pm);
%Pie chart
data = [Pfm,Protex,Pcomb,PexhExerg];
labels = {'Mechanical losses ( MW)','Mechanical exergetic power ( MW)','Primary exergy input ( MW)','Exhaust exergy losses ( MW)'};
pie(data,labels);

%Functions Masse Flow Rate
function MFair = massFlowRateAir(x,Pe,k_meca,lambda,pouvoirComburivore,enthalpieVector)
MFair=Pe*10^(3) - x*(1+(1/(lambda*pouvoirComburivore)))*(enthalpieVector(3)-enthalpieVector(4))*(1-k_meca)+x*(enthalpieVector(2)-enthalpieVector(1))*(1+k_meca);
end

function MFcomb = massFlowCombustible(pouvoirComburivore,lambda,MFair)
MFcomb=MFair/(lambda*pouvoirComburivore);
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
function Pexh = puissanceExhaust(lambda,pouvoirComburivore,enthalpieVector,MFair,MFfumee)
Pexh=(MFfumee*(1+1/(lambda*pouvoirComburivore))*enthalpieVector(4)-MFair*enthalpieVector(1))*10^(-3);
end

function Pm = puissanceMotrice(enthalpieVector,MFair,MFfumee)
Pm=(MFfumee*(enthalpieVector(3)-enthalpieVector(4)) - MFair*(enthalpieVector(2)-enthalpieVector(1)))*10^(-3);
end

function Pfm = puissanceFrottementMecanique(k_meca,enthalpieVector,MFair,MFfumee)
Pfm=k_meca*(MFfumee*(enthalpieVector(3)-enthalpieVector(4)) + MFair*(enthalpieVector(2)-enthalpieVector(1)))*10^(-3);
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

function PexhExerg = puissanceExhaustExerg(lambda,pouvoirComburivore,exergieVector,MFair,MFfumee)
PexhExerg =(MFfumee*(1+1/(lambda*pouvoirComburivore))*exergieVector(4)-MFair*exergieVector(1))*10^(-3);
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

function eta_cyclEn = rdtCyclEnergetique(lambda,pouvoirComburivore,enthalpieVector)
eta_cyclEn=((1+1/(lambda*pouvoirComburivore))*(enthalpieVector(3)-enthalpieVector(4))- (enthalpieVector(2)-enthalpieVector(1)))/((1+1/(lambda*pouvoirComburivore))*enthalpieVector(3)-enthalpieVector(2));
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

function eta_meca = rdtMeca (Pfm,Pm)
eta_meca=1-Pfm/Pm
end