%{
%%%%%%%%%%%%%%%%%%
Ref: Xin J, Gan L, Jiao L, Lai C. Accurate Density Calculation for Molten Slags in SiO2–Al2O3–CaO–MgO Systems. ISIJ International. 2017:ISIJINT-2017.

Notes: Only for SiO2-Al2O3-CaO-MgO

Output in g/cm^3
%%%%%%%%%%%%%%%%%%
%}


function [Rho_Xin] =  Rho_Xin(T,XSiO2,XTiO2,XAl2O3,XCr2O3,XFeO,XMgO,XMnO,XCaO,XNa2O,XK2O,XLi2O,XCaF2,XZrO2,XB2O3,XCrO,XNiO,XFe2O3,XBaO,XSrO)

V_ex = (5.64e-6)*XSiO2*XAl2O3...
    +(0.80e-6)*XSiO2*XCaO...
    +(2.00e-6)*XSiO2*XMgO...
    +(18.45e-6)*XAl2O3*XCaO...
    +(6.53e-6)*XAl2O3*XMgO...
    -(2.30e-6)*XCaO*XMgO;

V_m = XSiO2*V_SiO2(T)...
    +XAl2O3*V_Al2O3(T)...
    +XCaO*V_CaO(T)...
    +XMgO*V_MgO(T)...
    +V_ex;

MoleFraction =  [XSiO2, XTiO2, XAl2O3, XCr2O3,   XFeO,   XMgO, XMnO,   XCaO, XNa2O,  XK2O, XLi2O, XCaF2,  XZrO2, XB2O3, XCrO,  XNiO, XFe2O3,   XBaO,   XSrO]';
MolWeight = [  60, 79.9, 101.9,   152, 71.85,  40.3,  71, 56.08,   62, 94.2, 29.8,   78, 123.2, 69.8,  68, 74.7, 159.7, 153.3, 103.6]';
TotalMass = sum(MolWeight.*MoleFraction);
Rho_Xin = TotalMass/(V_m*1e6);

%% Functions for Temp Corrected Volume
    function V_SiO2 = V_SiO2(T)
        %Vi,R + DVi/DT(T-T_R)
        V_SiO2 = 26.312e-6 +0.740e-9 *(T-1773);
    end

    function V_Al2O3 = V_Al2O3(T)
        V_Al2O3 = 28.7e-6 + 10.108e-9 *(T-1773);
    end

    function V_CaO = V_CaO(T)
        V_CaO = 18.031e-6 + 1.014e-9 *(T-1773);
    end

    function V_MgO = V_MgO(T)
        V_MgO = 12.076e-6 + 0.683e-9 *(T-1773);
    end
end
