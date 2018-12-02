
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