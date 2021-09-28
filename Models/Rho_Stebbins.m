%{
Ref: Mills, K. C. (2011). The Estimation Of Slag Properties. Southern African Pyrometallurgy 2011 International Conference, March, 1–52. http://www.pyrometallurgy.co.za/KenMills/KenMills.pdf
 
DV_DT needs a better reference 

Output in g/cm^3
%}

function [Rho_Stebbins] =  Rho_Stebbins(T,XSiO2,XTiO2,XAl2O3,XCr2O3,XFeO,XMgO,XMnO,XCaO,XNa2O,XK2O,XLi2O,XCaF2,XZrO2,XB2O3,XCrO,XNiO,XFe2O3,XBaO,XSrO)

MoleFraction =  [XSiO2, XTiO2, XAl2O3, XCr2O3,   XFeO,   XMgO, XMnO,   XCaO, XNa2O,  XK2O, XLi2O, XCaF2,  XZrO2, XB2O3, XCrO,  XNiO, XFe2O3,   XBaO,   XSrO]';
MolWeight = [  60, 79.9, 101.9,   152, 71.85,  40.3,  71, 56.08,   62, 94.2, 29.8,   78, 123.2, 69.8,  68, 74.7, 159.7, 153.3, 103.6]';
feed_mol_frac_steb = MoleFraction(1:10)./(sum(MoleFraction)-MoleFraction(4)-MoleFraction(7)-sum(MoleFraction(11:end))); %renormalize because Cr and Mn aren't accounted for in the model
    
% Model
        % %         SiO2,  TiO2, Al2O3, Cr2O3,   FeO,    MgO, MnO,    CaO,   Na2O,    K2O
            a  =  [26.67, 24.42, 37.89,     0, 14.50,  12.39,   0,  17.68,  30.54,  48.92 ]'; % Don't forget to transponse!
            b  =  [-0.42,  0.81,  1.02,     0,  4.16,   2.22,   0,   3.97,   7.43,  11.99 ]'/1000;

    Rho_Stebbins = sum(feed_mol_frac_steb.*MolWeight(1:10))./(sum(feed_mol_frac_steb.*a)+(T-1873)*sum(feed_mol_frac_steb.*b)); %g/cm^3
    
end
