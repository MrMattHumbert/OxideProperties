%{
%%%%%%%%%%%%%%%%%%
Ref: Mills, K. C. (2011). The Estimation Of Slag Properties. Southern African Pyrometallurgy 2011 International Conference, March, 1–52. http://www.pyrometallurgy.co.za/KenMills/KenMills.pdf 

Notes: Reproduced from Mills

Output in dPa-s
%%%%%%%%%%%%%%%%%%
%}

function [ViscosityUrbain] =  Vis_Urbain(T,XSiO2,XTiO2,XAl2O3,XCr2O3,XFeO,XMgO,XMnO,XCaO,XNa2O,XK2O,XLi2O,XCaF2,XZrO2,XB2O3,XCrO,XNiO,XFe2O3,XBaO,XSrO)

    NormUrb = 1+0.5*(XFe2O3+XCr2O3)+XTiO2+XZrO2+2*XCaF2;
    XgUrb = XSiO2/NormUrb;
    XmUrb = (2*(XTiO2+XZrO2)+XFeO+XMgO+XMnO+XCaO+XNa2O+XK2O+XLi2O+3*XCaF2+XCrO+XNiO+XBaO+XSrO+3*0.6*(XCr2O3+XFe2O3))/NormUrb;
    XaUrb = (XAl2O3+XB2O3+0.4*(XCr2O3+XFe2O3))/NormUrb;
    alphaUrb = XmUrb/(XmUrb+XaUrb);
    B0 = 13.8+39.9355*alphaUrb-44.049*alphaUrb^2;
    B1 = 30.481-117.15*alphaUrb+129.9978*(alphaUrb)^2;
    B2 = -40.9429+234.0486*alphaUrb-300.04*(alphaUrb)^2;
    B3 = 60.7619-153.9276*alphaUrb+211.1616*(alphaUrb)^2;
    B = B0+B1*XgUrb+B2*XgUrb^2+B3*XgUrb^3;
    
    Ca_mod =(13.2+41.5*alphaUrb-45*(alphaUrb)^2)+(30.5-117.2*alphaUrb+130*(alphaUrb)^2)*XgUrb+(-40.4+232.1*alphaUrb-298.6*(alphaUrb)^2)*XgUrb^2+(60.8-156.4*alphaUrb+213.6*(alphaUrb)^2)*XgUrb^3;
    Mg_mod =(13.2+15.9*alphaUrb-18.6*(alphaUrb)^2)+(30.5-54.1*alphaUrb+33*(alphaUrb)^2)*XgUrb+(-40.4+138*alphaUrb-112*(alphaUrb)^2)*XgUrb^2+(60.8-99.8*alphaUrb+97.6*(alphaUrb)^2)*XgUrb^3;
    Mn_mod =(13.2+20*alphaUrb-25.6*(alphaUrb)^2)+(30.5+26*(alphaUrb)-56*(alphaUrb)^2)*XgUrb+(-40.4-110.3*alphaUrb+186.2*(alphaUrb)^2)*XgUrb^2+(60.8+64.3*alphaUrb-104.6*(alphaUrb)^2)*XgUrb^3;
    XmnUrb = XFeO+XMnO+XCrO+XCr2O3+XNiO+XFe2O3;
    XcUrb = XCaO+XBaO+XSrO+XNa2O+XK2O+XLi2O;
    B_global =(Ca_mod*XcUrb+Mg_mod*XMgO+Mn_mod*XmnUrb)/(XcUrb+XmnUrb+XMgO);
    AUrb = exp(-(0.29*B_global+11.57));
    ViscosityUrbain = AUrb*T*exp(1000*B/T);


end