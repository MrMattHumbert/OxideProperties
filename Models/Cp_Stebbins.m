%{
Ref: Ref: Mills, K. C. (2011). The Estimation Of Slag Properties. Southern African Pyrometallurgy 2011 International Conference, March, 1–52. http://www.pyrometallurgy.co.za/KenMills/KenMills.pdf 
Notes: This data is dated. 

Output in J/(K-gfw) gfw==mol
%}

function [Cp_T_Steb, Del_H_T_Steb, Cp_T_L_Steb] = Cp_Stebbins(T,XSiO2,XTiO2,XAl2O3,XCr2O3,XFeO,XMgO,XMnO,XCaO,XNa2O,XK2O,XLi2O,XCaF2,XZrO2,XB2O3,XCrO,XNiO,XFe2O3,XBaO,XSrO)

feed_mol = [XSiO2,XTiO2,XAl2O3,XCr2O3,  XFeO,  XMgO,XMnO,  XCaO,XNa2O, XK2O,XLi2O,XCaF2, XZrO2,XB2O3,XCrO,XNiO, XFe2O3,  XBaO,  XSrO]';
    

    a_steb = [ 66.345, 33.851, 91.404,    0, 40.949, 32.244,  0, 46.677, 69.067, 107.194, ]';
    b_steb = [ 0.7797, 6.4572, 4.4940,    0, 2.9349, 2.7288,  0, 0.3565, 1.8603, -3.2194, ]'/100;
    c_steb = [-28.003,  4.470,-21.465,    0,-7.6986, 1.7549,  0,-1.9322, 2.9101, -28.929, ]'/1e-5;
    
    feed_mol_frac_steb = feed_mol(1:10)./(sum(feed_mol)-feed_mol(4)-feed_mol(7)-sum(feed_mol(11:end))); %renormalize because Cr and Mn aren't accounted for in the model
    
    Cp_T_Steb = (sum(a_steb.*feed_mol_frac_steb)+ sum(b_steb.*feed_mol_frac_steb)*T + sum(c_steb.*feed_mol_frac_steb)*T^(-2)); %g/K-mol
    Del_H_T_Steb = sum(a_steb.*feed_mol_frac_steb)*(T-298)+ 0.5*sum(b_steb.*feed_mol_frac_steb)*T^2 -0.5*sum(b_steb.*feed_mol_frac_steb)*(298)^2 + sum(c_steb.*feed_mol_frac_steb)/T - sum(c_steb.*feed_mol_frac_steb)/298;
    
% Stebbins J/k-gfw gfw==mol LIQUID
% Data
%             SiO2,   TiO2,  Al2O3, Cr2O3,    FeO,    MgO,   MnO,    CaO,   Na2O,     K2O,
    Cp_l_Steb = [ 80.0, 111.8, 157.6, 0, 78.9,99.7,0,99.9,102.3,97 ]';
    Cp_T_L_Steb = sum(feed_mol_frac_steb.*Cp_l_Steb);
    
end