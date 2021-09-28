%{
%%%%%%%%%%%%%%%%%%
Ref: Mills, K. C. (2011). The Estimation Of Slag Properties. Southern African Pyrometallurgy 2011 International Conference, March, 1–52. http://www.pyrometallurgy.co.za/KenMills/KenMills.pdf 
Ref: Kekkonen, M., Oghbasilasie, H., & Louhenkilpi, S. (2012). Viscosity models for molten slags. Aalto University Publication Series SCIENCE + TECHNOLOGY, 12, 38.

Notes: Mills worksheet has errors for the Riboud model Use Kekkonen

Output in dPa-s
%%%%%%%%%%%%%%%%%%
%}

function [ViscosityRiboud] =  Vis_Riboud(T,XSiO2,XTiO2,XAl2O3,XCr2O3,XFeO,XMgO,XMnO,XCaO,XNa2O,XK2O,XLi2O,XCaF2,XZrO2,XB2O3,XCrO,XNiO,XFe2O3,XBaO,XSrO)

    XsRib = XSiO2+XTiO2+XZrO2;
    
    XalkRib = XNa2O+XK2O+XLi2O+XBaO+XSrO;
    XaRib = XAl2O3+XB2O3;
    XnbRib = XCaO+XMgO+XFeO+XMnO+XCr2O3+XCrO+XNiO+XFe2O3;
    A_Rib = exp(-19.81+1.73*XnbRib+5.82*XCaF2+7.02*XalkRib-35.76*XaRib);
    B_Rib = 31140-23896*XnbRib-46356*XCaF2-39159*XalkRib+68833*XaRib;
    ViscosityRiboud = A_Rib*T*exp(B_Rib/T);

end