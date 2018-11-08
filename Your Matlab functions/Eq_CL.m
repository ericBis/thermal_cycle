function F = Eq_CL(x,Phi_out,ha_in,he2,he3,he4,xa_in,Pw)
% Inputs
% x: same form that the ouput F
% Phi_out: [-] Maximum relative humidity of air at the cooling 
%                         tower outlet.
% ha_in: [J/kg] air enthalpy at the tower inlet
% he2: [J/kg] water enthalpy at point 2
% he3: [J/kg] water enthalpy at point 3
% he4: [J/kg] water enthalpy at point 4
% xa_in: [\] absolute humidity of the air at the tower inlet
% P_w: [J/kg] Heat power to dissipate


% Outputs
% F(1)= ha_out
% F(2)= xa_out
% F(3)= m_as
% F(4)= me1
% F(5)= me2=me3
% F(6)= me4
% F(7)= he1

F=zeros(7,0);
[~,~,~,ha_out,~,~,~] = Psychrometrics('w',x(2),'phi',Phi_out*100);

F(1)=x(1)-ha_out;
F(2)= x(3)*(ha_in-x(1)) + x(5)*(he3-x(7)) + x(6)*x(7);
F(3)= x(3)*(x(2)-xa_in) - x(6);
F(4)= Pw - x(5)*he3 + x(4)*x(7);
F(5)= x(4)+ x(6) - x(5);
F(6)= Pw - x(5)*he3 + x(4)*x(7);
F(7)= x(4)*x(7)+ x(6)*he4 - x(5)*he2;
end

