%{
%%%%%%%%%%%%%%%%%%
Ref: Giordano, D., & Dingwell, D. B. (2003). Non-Arrhenian multicomponent melt viscosity: A model. Earth and Planetary Science Letters, 208(3–4). https://doi.org/10.1016/S0012-821X(03)00042-6
Ref: Erratum for the paper is very important!

Notes: NBO/T in the presented model was originally from Mysen (B.O. Mysen, Structure
and Properties of Silicate Melts, Elsevier, Amsterdam 1988, 354 pp.),
Herein the values are taken from Mills

Output in dPa-s
%%%%%%%%%%%%%%%%%%
%}

function [Eta_NBO,Eta_SM] =  Vis_Giordano(T,XSiO2,XTiO2,XAl2O3,XCr2O3,XFeO,XMgO,XMnO,XCaO,XNa2O,XK2O,XLi2O,XCaF2,XZrO2,XB2O3,XCrO,XNiO,XFe2O3,XBaO,XSrO)

% Supporting Functions
[~,~,~,NBO_T2] = Q(T,XSiO2,XTiO2,XAl2O3,XCr2O3,XFeO,XMgO,XMnO,XCaO,XNa2O,XK2O,XLi2O,XCaF2,XZrO2,XB2O3,XCrO,XNiO,XFe2O3,XBaO,XSrO);

a1 = -0.15139 - 1129.19/(T-273) - 1381914/(T-273)^2 + (1.29e9)/(T-273)^3 ;
a2 = -0.00071 - 3.47074/(T-273) + 5720.781/(T-273)^2 - 2061030/(T-273)^3 ;
a3 = -5.44516 + 9309.88/(T-273) - 3390935/(T-273)^2 + 3144300695/(T-273)^3; 

Eta_NBO = (10^(a1*log(NBO_T2-a2)+a3))*10 ;

c1 = (-17.80106 + 0.01808103*(T-273))/(1-(2.2869e-3)*(T-273));
c2 = 1/(0.02532 + 2.5124*exp((-6.3679e-3)*(T-273))+40.4562e-6*(T-273));
c3 = (1-1.6569e-3*(T-273))/(0.017954-63.90597e-6*(T-273));
SM = (XNa2O + XK2O + XCaO + XMgO + XMnO + (XFeO + XFe2O3)/2)*100; 

Eta_SM = (10^(c1 + c2*c3/(c3 + SM)))*10;


end