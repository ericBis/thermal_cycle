
tspanc=[288.15,702.15];
tspancomb=[702.15,1673.15];
tspand=[1673.15,927.15];
lambda=2.3;
s1=0.054;
X=0;
Y=4;
[tc,sc]=ode45(@(tc,sc) (1-0.9)*CpAir(tc)/tc,tspanc,s1);
scc=sc(end);
[tcomb,scomb]=ode45(@(tcomb,scomb) dsdt(tcomb,Y,X,lambda),tspancomb,scc);
[td,sdd]=ode45(@(td,sd) (-(1-0.9)/0.9)*CpFumee(td,Y,X,lambda)/td,tspand,scomb(end));
abscisseEntrop=[sc',scomb',sdd'];
ordonneeTemperature=[tc',tcomb',td'];
plot(abscisseEntrop,ordonneeTemperature);
% vecteurentrop=[0.0540,0.1418,1.3375,1.4181]
% sc(1)=0.054
% sc(end)=0.1461
% scomb(end)=1.2231
% sdd(end)=1.3065 valeur differentes car formule differentes + pas de
% temps de ODE45
%A CORRIGER ET FAIRE LE H-S

function DSDT= dsdt(t,Y,X,lambda)
if t<=702.3
DSDT=CpAir(t)/t;
else
    DSDT=CpFumee(t,Y,X,lambda)/t;
end
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

function RstarAir = constanteRair()
R=8.314e-3;%[kJ/mol.K]
MmO2=32e-3;%[kg/mol]
MmN2=28e-3;%[kg/mol]
propO2=0.21;%Propotion en mole
propN2=0.79;%Propotion en mole
RstarAir=R/(MmO2*propO2+MmN2*propN2);%[kJ/kgAirK]
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
