%%%%%%%%%%%%%%%%%%APPEL FONCTION ETAT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fun1 = @calculNaT2;
x00 = [1.2,295];
solutionNaetT2 = lsqnonlin(fun1,x00); %Determination de na et T2
%%%
fun00 = @calculNgT4;
x0 =[700, 1.1];%x0(1)=T4,x0(2)=lambda
X=0;
Y=4;
solutionT4etLambda = fsolve(@(x)fun00(x,Y,X),x0) %Determination de T4 et lambda
%%%
T1=288;
T3=1673;
TemperatureEtats=[T1,solutionNaetT2(2),T3,solutionT4etLambda(1)]
enth = enthalpy(TemperatureEtats)%Determination du vect enthalpie
%%%
pressions = [1,18,17.1,1];
entrop = entropy(TemperatureEtats,pressions)
%%%
exergie = exergy(enth,entrop)
%%%
matriceEtats=ones(4,4);
matriceEtats(1,:)=TemperatureEtats;
matriceEtats(2,:)=pressions;
matriceEtats(3,:)=enth;
matriceEtats(4,:)=entrop;
matriceEtats(5,:)=exergie;
etat=matriceEtats
%%%
LHV=PouvoirCalorifiqueInf(Y,X);
ec=exergyComb;
lambda=solutionT4etLambda(2);
cpFumeeAt400K=cpMoyenFumee(400,927,Y,X,lambda)/(solutionT4etLambda(1)-400); %kJ/kgK
%%%
% plot()
%%%%%%%%%%%%%%%%%%%%%%%%%%%PARTIE ETAT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function enth = enthalpy(T)
%TemperatureEtats(vecteur LIGNE avec T1,T2,T3,T4) EN KELVIN
%cpMoyen est un vecteur de cpMoyen calculer correctement pour chaque etat
%et calculer à chaque fois par rapport à la référence T0
%(cpMoyenAir1,cpMoyenAir2,cpMoyenFumee3,cpMoyenFumee4) LIGNE
T0=273;%EN KELVIN
Y=4;
X=0;
lambda=2.3;
cpmoyenAirEtFumee=ones(1,4);
cpmoyenAirEtFumee(1)=cpMoyenAir(T0,T(1));
cpmoyenAirEtFumee(2)=cpMoyenAir(T0,T(2));
cpmoyenAirEtFumee(3)=cpMoyenFumee(T0,T(3),Y,X,lambda);
cpmoyenAirEtFumee(4)=cpMoyenFumee(T0,T(4),Y,X,lambda);
enth=cpmoyenAirEtFumee.*(T-T0);%rendre vecteur de meme taille que "TemperatureEtats"
end

function entrop = entropy(Temp,pressions)
entrop=ones(1,4);
p0=1;
T0=273;
X=0;
Y=4;
lambda=2.3;

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
T0=288;
exergie=ones(1,4);
exergie(1)=(h(1)-h0)-T0*(s(1)-s0);
exergie(2)=(h(2)-h0)-T0*(s(2)-s0);
exergie(3)=(h(3)-h0)-T0*(s(3)-s0);
exergie(4)=(h(4)-h0)-T0*(s(4)-s0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PARTIE COMBUSTION%%%%%%%%%%%%%%%%%%%%%%
function eqq = calculNgT4(x,Y,X)
eta_PiT=0.9;%x(1)=T4 x(2)=lambda
k_cc=0.95;
Tcomb=298;
r=18;
T0=273;
T2=702;
T3=1673;
eqq=ones(2,1);
cpMoyenComb=(35.26)/(MasseMolaireCombustible(Y,X));
eqq(1,1)=log(x(1)/T3) - (eta_PiT*constanteRfumee(Y,X,x(2))/cpMoyenFumee(x(1),T3,Y,X,x(2)))*log(1/(k_cc*r));
eqq(2,1)=(T3-T0)-(PouvoirCalorifiqueInf(Y,X)+x(2)*pouvoirComburivore(Y,X)*cpMoyenAir(T0,T2)*(T2-T0)+cpMoyenComb*(Tcomb-T0))/(((pouvoirComburivore(Y,X))*x(2)+1)*cpMoyenFumee(T0,T3,Y,X,x(2)));
end

function cpMFumee = cpMoyenFumee(Ta,Tb,Y,X,lambda)%Valable pour les etats 2,3,4 et tient compte de la compostition de l'air
fun=@CpFumee;
cpMFumee=(integral(@(T)fun(T,Y,X,lambda),Ta,Tb))/(Tb-Ta);
end

function cpF = CpFumee(a,y,x,lambda)%Valable pour les etats 3 et 4, utiliser janaf directement. Compositions de l'air:O2 N2, CO2, H2O
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

propO2=(lambda*(1+y/4-x/2)-1+x/2-y/4)/(4.76*lambda*(1+y/4-x/2)+x/2+y/4);%Propotion molO2/molFumee
propCO2=1/(4.76*lambda*(1+y/4-x/2)+x/2+y/4);%Propotion molCO2/molFumee
propH2O=(y/2)/(4.76*lambda*(1+y/4-x/2)+x/2+y/4);%Propotion molH2O/molFumee
propN2=(3.76*(lambda*(1+y/4-x/2)))/(4.76*lambda*(1+y/4-x/2)+x/2+y/4);%Propotion molN2/molFumee

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

function pci =PouvoirCalorifiqueInf(Y,X)%tiree d'une formule empirique
formuleEmpirique=393400+102250*Y-(X/(1+0.5*Y))*(111000+102250*Y);
pci=formuleEmpirique/(MasseMolaireCombustible(Y,X));%donne le PCI en kJ/kgComb de type CHyOx c'est un scalaire
end

function pcs = PouvoirCalorifiqueSup()
pcs = 55695;%donne le PCS en kJ/kgComb
end

function ec = exergyComb()
ec = 52215;%donne ec du CH4 en kJ/kgComb
end

function ma1 = pouvoirComburivore(y,x)%Sachant que x et y sont les proprietes du fuel
ma1=(1+(y-2*x)/4)*(32+3.76*28.15)/(12+y+16*x);%[kgair/kgcomb] c'est un scalaire
end

function MMComb = MasseMolaireCombustible(y,x)%Pour un combustible de type CHyOx
MmO=16;%[kg/kmolO]
MmH=1;%[kg/kmolH]
MmC=12;%[kg/kmolC]
MMComb=MmC+y*MmH+x*MmO;%[kgComb/kmolComb] c'est un scalaire
end

%%%%%%%%%%%%%%%%%%%%%%%%PARTIE ADMISSION AIR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function eq= calculNaT2(x)
eta_PiC=0.9;%x(1)=na et x(2)=T2
T1=288;% TOUTE LES TEMP DU CODE SONT EN KELVIN
r=18;
eq=ones(2,length(x));
eq(1,1)= (x(1)-1)/x(1) - (1/eta_PiC)*(constanteRair/cpMoyenAir(T1,x(2)));
eq(2,1)= x(2)/T1 -r^((x(1)-1)/x(1));
end

function RstarAir = constanteRair()
R=8.314e-3;%[kJ/mol.K]
MmO2=32e-3;%[kg/mol]
MmN2=28e-3;%[kg/mol]
propO2=0.21;%Propotion en mole
propN2=0.79;%Propotion en mole
RstarAir=R/(MmO2*propO2+MmN2*propN2);%[kJ/kgAirK]
end

function cpMAir = cpMoyenAir(Ta,Tb)%Valable pour les etats 1 et 2 et tient compte de la compostition de l'air
fun=@CpAir;
cpMAir=(integral(fun,Ta,Tb))/(Tb-Ta);%[kJ/kgAirK]
end

function cpA = CpAir(a)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%JANAF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out=janaf(prop, spec,T)
% ---------------------------------------------------------------------
% function out=janaf(prop, spec, T)                 |   Version 1.01
% ---------------------------------------------------------------------
% Calculates JANAF curve fit according to JANAF virial equation
% Output is calculated in SI-units
%
% prop = 'c' for standard state specific heat
%      = 'h' for standard state enthalpy
%      = 's' for standard state entropy
% spec = 'CO2', 'H2O', 'CO', 'H2', 'O2', 'N2'
% T    = Temperature, vector allowed
% 
% JANAF.mat required (contains the coefficients)
% ---------------------------------------------------------------------
% Last Change: 2003-07-18           |   (c)2003, Stefan Billig, Delphi

% check for correct syntax
if nargin~=3
    help janaf
    % end function
    return
end

% load coefficient table
CO = [2.98410000000000,0.00148910000000000,-5.79000000000000e-07,1.03650000000000e-10,-6.93540000000000e-15,-14245,6.34790000000000;3.71010000000000,-0.00161910000000000,3.69240000000000e-06,-2.03200000000000e-09,2.39530000000000e-13,-14356,2.95550000000000;];
CO2 = [4.46080000000000,0.00309820000000000,-1.23930000000000e-06,2.27410000000000e-10,-1.55260000000000e-14,-48961,-0.986300000000000;2.40080000000000,0.00873510000000000,-6.60710000000000e-06,2.00220000000000e-09,6.32740000000000e-16,-48378,9.69510000000000;];
H2 = [3.10020000000000,0.000511190000000000,5.26440000000000e-08,-3.49100000000000e-11,3.69450000000000e-15,-877.380000000000,-1.96290000000000;3.05740000000000,0.00267650000000000,-5.80990000000000e-06,5.52100000000000e-09,-1.81230000000000e-12,-988.900000000000,-2.29970000000000;];
H2O = [2.71680000000000,0.00294510000000000,-8.02240000000000e-07,1.02270000000000e-10,-4.84720000000000e-15,-29906,6.63060000000000;4.07010000000000,-0.00110840000000000,4.15210000000000e-06,-2.96370000000000e-09,8.07020000000000e-13,-30280,-0.322700000000000;];
N2 = [2.89630000000000,0.00151550000000000,-5.72350000000000e-07,9.98070000000000e-11,-6.52240000000000e-15,-905.860000000000,6.16150000000000;3.67480000000000,-0.00120820000000000,2.32400000000000e-06,-6.32180000000000e-10,-2.25770000000000e-13,-1061.20000000000,2.35800000000000;];
O2 = [3.62200000000000,0.000736180000000000,-1.96520000000000e-07,3.62020000000000e-11,-2.89460000000000e-15,-1202,3.61510000000000;3.62560000000000,-0.00187820000000000,7.05550000000000e-06,-6.76350000000000e-09,2.15560000000000e-12,-1047.50000000000,4.30530000000000;];
MolWeight.O2 = 31.998000000000000;
MolWeight.CO2 = 44.008000000000000;
MolWeight.H2 = 2.015940000000000;
MolWeight.H2O = 18.014940000000000;
MolWeight.N2 = 28.014000000000000;
MolWeight.CO = 28.0090000000000;

z=1;
% determine molecular weight from table
MWeight=eval(['MolWeight.' spec]);

for i=1:length(T)
    % choose temperature range vector
    if (T(i)>1000 & T(i)<=5000)
        eval(['ai=' spec '(1,:);'])
        out(z)=calc(ai, T(i), prop, MWeight);
        z=z+1;
    elseif (T(i)>=300 & T(i)<=1000)
        eval(['ai=' spec '(2,:);'])
        out(z)=calc(ai, T(i), prop, MWeight);
        z=z+1;
    else
        sprintf(['Temperature ' num2str(T(i)) 'K not between 300K and 5000K!'])
    end
end
end

% out = zeros(size(T));
% ind1 = T>1000 & T<=5000;
% eval(['ai=' spec '(1,:);']);
% out(ind1) = calc(ai, T(ind1), prop, MWeight);
% 
% ind2 = T>=300 & T<=1000;
% eval(['ai=' spec '(2,:);']);
% out(ind2) = calc(ai, T(ind2), prop, MWeight);
% 
% ind3 = T>5000 & T<300;
% if ~isempty(T(ind3))
%     sprintf(['Temperatures ' num2str(T(ind3)) 'K not between 300K and 5000K!']);
% end

%----------------------------------------------------------------------
function out=calc(ai, T, prop, MWeight)

R=8.314472;
% calculate standard state value
switch prop
    case 'c'
        out=(ai(1)+ai(2)*T+ai(3)*T.^2+ai(4)*T.^3+ai(5)*T.^4)*R/MWeight;
    case 'h'
        out=(ai(1)+ai(2)/2*T+ai(3)/3*T.^2+ai(4)/4*T.^3+ai(5)/5*T.^4+ai(6)/T).*T*R/MWeight;
    case 's'
        out=(ai(1)*log(T)+ai(2)*T+ai(3)/2*T.^2+ai(4)/3*T.^3+ai(5)/4*T.^4+ai(7))*R/MWeight;
end
end
