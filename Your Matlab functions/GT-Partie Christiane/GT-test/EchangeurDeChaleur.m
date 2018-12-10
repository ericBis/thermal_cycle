%%%%Ajout de l'echangeur de chaleur%%%%
% Temperatures conservees: t1,t2,t3,t4
% Temperatures recherchee: t2R,t5
NTU=[0,1,2,4,8];
coef=4./(1+3*NTU)
options=struct;
options.eta_PiC=0.9;
options.eta_PiT=0.9;
options.k_cc=0.95;
options.k_mec=0.015;
options.T_3=1400;
options.T_0=15;
T3=273.15+options.T_3;
T1=273.15+options.T_0;
t4=654;
T4=273.15+t4;
t2=429;
T2=273.15+t2;
lambda=2.3;
X=0;
Y=4;
%%%%Calcul de T2R en fonction du NTU%%%%
T2R=(T4.*NTU+T2)./(1+NTU)
t2R=T2R-273.15
% %%%%Calcul de T5 en fonction d'un bilan d'energie%%%%
% fun000=@calculTemperatureT5;
% x000=85;%[K]
% fsolve
% T5 = fsolve(@(x)fun000(x),x000);
% function calculT5 = calculTemperatureT5(x)
% calculT5=MFair*cpMoyenAir(T2,T2R)*(T2R-T2)-MFfumee*cpMoyenFumee(T4,x,Y,X,lambda)*(T4-x);
% end
