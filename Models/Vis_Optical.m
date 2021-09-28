%{
Ref: Mills, K. C. (2011). The Estimation Of Slag Properties. Southern African Pyrometallurgy 2011 International Conference, March, 1–52. http://www.pyrometallurgy.co.za/KenMills/KenMills.pdf 

Notes: Reproduced from Mills

Output in dPa-s
%}

function [ViscosityOptical] =  Vis_Optical(T,XSiO2,XTiO2,XAl2O3,XCr2O3,XFeO,XMgO,XMnO,XCaO,XNa2O,XK2O,XLi2O,XCaF2,XZrO2,XB2O3,XCrO,XNiO,XFe2O3,XBaO,XSrO)

% Supporting Functions
Lam_Corr = OpticalBasicity(T,XSiO2,XTiO2,XAl2O3,XCr2O3,XFeO,XMgO,XMnO,XCaO,XNa2O,XK2O,XLi2O,XCaF2,XZrO2,XB2O3,XCrO,XNiO,XFe2O3,XBaO,XSrO);


% Model
A_opt = exp(-232.69*Lam_Corr^2 +357.32*Lam_Corr -144.17);
B_opt = 1000*exp(-1.77+2.88/Lam_Corr);
ViscosityOptical = 10*A_opt*exp(B_opt/T); 

end