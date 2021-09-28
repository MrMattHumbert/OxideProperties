%{
Ref: Mills, K. C. (2011). The Estimation Of Slag Properties. Southern African Pyrometallurgy 2011 International Conference, March, 1–52. http://www.pyrometallurgy.co.za/KenMills/KenMills.pdf
Ref: Mills, K. C., Karagadde, S., Lee, P. D., Yuan, L., & Shahbazian, F. (2016). Calculation of physical properties for use in models of continuous casting process-Part 1: Mould slags. ISIJ International, 56(2). https://doi.org/10.2355/isijinternational.ISIJINT-2015-364

Notes: 

Output in J/(K-mol)
%}

function [Cp_Mills] =  Cp_Mills(T,XSiO2,XTiO2,XAl2O3,XCr2O3,XFeO,XMgO,XMnO,XCaO,XNa2O,XK2O,XLi2O,XCaF2,XZrO2,XB2O3,XCrO,XNiO,XFe2O3,XBaO,XSrO)

feed_mol_frac = [XSiO2,XTiO2,XAl2O3,XCr2O3,  XFeO,  XMgO,XMnO,  XCaO,XNa2O, XK2O,XLi2O,XCaF2, XZrO2,XB2O3,XCrO,XNiO, XFe2O3,  XBaO,  XSrO]';

          %    SiO2,    TiO2, Al2O3, Cr2O3,    FeO,    MgO,    MnO,    CaO,  Na2O,    K2O,  Li2O,CaF2,     ZrO2,    B2O3,   CrO,   NiO, Fe2O3,  BaO,  SrO
a_mills = [86.98536,111.7128,146.44,146.44,76.5672,90.3744,79.9144,80.7512,92.048,74.0568,96.232,96.232,112.968,126.7752,75.312,75.312,146.44,83.68,83.68]';
Cp_Mills = sum(a_mills.*feed_mol_frac);  

end