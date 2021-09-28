%{
Ref: Mills, K. C. (2011). The Estimation Of Slag Properties. Southern African Pyrometallurgy 2011 International Conference, March, 1–52. http://www.pyrometallurgy.co.za/KenMills/KenMills.pdf
Ref: Mills, K. C., & Keene, B. J. (1987). Physical properties of BOS slags. International Materials Reviews, 32(1), 1–120. https://doi.org/10.1179/095066087790150296

Output in g/cm^3
%}

function [Rho_Keen] =  Rho_Keen(T,XSiO2,XTiO2,XAl2O3,XCr2O3,XFeO,XMgO,XMnO,XCaO,XNa2O,XK2O,XLi2O,XCaF2,XZrO2,XB2O3,XCrO,XNiO,XFe2O3,XBaO,XSrO)
%                SiO2, TiO2, Al2O3, Cr2O3,   FeO,   MgO, MnO,   CaO, Na2O,  K2O, Li2O, CaF2,  ZrO2, B2O3, CrO,  NiO, Fe2O3,   BaO,   SrO
MoleWeight =    [  60, 79.9, 101.9,   152, 71.85,  40.3,  71, 56.08,   62, 94.2, 29.8,   78, 123.2, 69.8,  68, 74.7, 159.7, 153.3, 103.6]';
MoleFraction = [XSiO2,XTiO2,XAl2O3,XCr2O3,  XFeO,  XMgO,XMnO,  XCaO,XNa2O, XK2O,XLi2O,XCaF2, XZrO2,XB2O3,XCrO,XNiO, XFe2O3,  XBaO,  XSrO]';
TotalMass = sum(MoleWeight.*MoleFraction); %Total mass of 1 mol

%Convert to mass fraction
mFeO = MoleWeight(5)*XFeO/TotalMass;
mMnO = MoleWeight(7)*XMnO/TotalMass;
mCr2O3 = MoleWeight(4)*XCr2O3/TotalMass;
mCrO = MoleWeight(15)*XCrO/TotalMass;
mNiO = MoleWeight(16)*XNiO/TotalMass;
mFe2O3 = MoleWeight(17)*XFe2O3/TotalMass;


%kg/m^3
Rho_Keen_kg = (2490 + 1200*(mFeO+mMnO+mCr2O3+mCrO+mNiO+mFe2O3))*(1-(T-1673)/10000); % 12 from Mills "Estimation of Slag Properties" 18 cited elsewhere

Rho_Keen = Rho_Keen_kg/1000;
    
end
