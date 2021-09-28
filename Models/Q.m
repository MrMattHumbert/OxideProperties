%{
%%%%%%%%%%%%%%%%%%
Ref: Mills, K. C. (2011). The Estimation Of Slag Properties. Southern African Pyrometallurgy 2011 International Conference, March, 1–52. http://www.pyrometallurgy.co.za/KenMills/KenMills.pdf

Notes: Does not incorportate CaF2 and SrO in the structural Model. 

Output is dimensionless
%%%%%%%%%%%%%%%%%%
%}

function [Q,Qf,NBO_T,NBO_T2] = Q(T,XSiO2,XTiO2,XAl2O3,XCr2O3,XFeO,XMgO,XMnO,XCaO,XNa2O,XK2O,XLi2O,XCaF2,XZrO2,XB2O3,XCrO,XNiO,XFe2O3,XBaO,XSrO)

XT = XSiO2+2*XAl2O3+0.45*2*XFe2O3+0.55*2*XCr2O3;
XNB_1 = 2*(XCaO+XMgO+XNa2O+XK2O+XLi2O+XFeO+XMnO+XCrO+XNiO+XBaO+2*(XTiO2+XZrO2)+3*XB2O3+0.45*4*XCr2O3+0.55*4*XFe2O3-XAl2O3-XCr2O3-XFe2O3);
NBO_T = XNB_1/XT;
Q = 4-NBO_T;
NBO_T2 = 2*(XCaO+XMgO+XNa2O+XK2O+XLi2O+XFeO+XMnO+XCrO+XNiO+XBaO+2*(XTiO2+XZrO2)+3*XB2O3+(0.73-0.113*Q)*4*XCr2O3+(0.67-0.086*Q)*4*XFe2O3-XAl2O3-XCr2O3-XFe2O3)/...
    (XSiO2+2*XAl2O3+(0.33+0.086*Q)*2*XFe2O3+(0.27+0.113*Q)*2*XCr2O3);
Qf = 4-NBO_T2;

end

