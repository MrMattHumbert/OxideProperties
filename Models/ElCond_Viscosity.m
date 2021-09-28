%{
Ref: Mills, K. C. (2011). The Estimation Of Slag Properties. Southern African Pyrometallurgy 2011 International Conference, March, 1–52. http://www.pyrometallurgy.co.za/KenMills/KenMills.pdf
 
Notes: Reproduced from Mills

Output in S/cm
%}

function [ElCond_Giordano, ElCond_Optical, ElCond_Riboud, ElCond_Urbain] = ElCond_Viscosity(T,XSiO2,XTiO2,XAl2O3,XCr2O3,XFeO,XMgO,XMnO,XCaO,XNa2O,XK2O,XLi2O,XCaF2,XZrO2,XB2O3,XCrO,XNiO,XFe2O3,XBaO,XSrO)

% Call supporting Functions
    [~,ViscosityGiordano] = Vis_Giordano(T,XSiO2,XTiO2,XAl2O3,XCr2O3,XFeO,XMgO,XMnO,XCaO,XNa2O,XK2O,XLi2O,XCaF2,XZrO2,XB2O3,XCrO,XNiO,XFe2O3,XBaO,XSrO);
    ViscosityOptical = Vis_Optical(T,XSiO2,XTiO2,XAl2O3,XCr2O3,XFeO,XMgO,XMnO,XCaO,XNa2O,XK2O,XLi2O,XCaF2,XZrO2,XB2O3,XCrO,XNiO,XFe2O3,XBaO,XSrO);
    ViscosityRiboud = Vis_Riboud(T,XSiO2,XTiO2,XAl2O3,XCr2O3,XFeO,XMgO,XMnO,XCaO,XNa2O,XK2O,XLi2O,XCaF2,XZrO2,XB2O3,XCrO,XNiO,XFe2O3,XBaO,XSrO);
    ViscosityUrbain = Vis_Urbain(T,XSiO2,XTiO2,XAl2O3,XCr2O3,XFeO,XMgO,XMnO,XCaO,XNa2O,XK2O,XLi2O,XCaF2,XZrO2,XB2O3,XCrO,XNiO,XFe2O3,XBaO,XSrO);

% Model

    if (XNa2O+XLi2O+XK2O)>0.005 % Method2
        Xna = 2*(XNa2O+XK2O+XLi2O);
        Xcat_tot = XCaO+XMgO+Xna+XFeO+XMnO+XCaF2+XCrO+XNiO+XBaO+XSrO+0.5*(XZrO2+XTiO2)+0.6667*(XCr2O3+XFe2O3);
        F = 0.15+(Xna/Xcat_tot*3.87);
        G = 1.1+1.77*Xna/Xcat_tot;
        
        ElCond_Giordano = exp((F-log(ViscosityGiordano))/G);
        ElCond_Optical  = exp((F-log(ViscosityOptical))/G);
        ElCond_Riboud = exp((F-log(ViscosityRiboud))/G);
        ElCond_Urbain = exp((F-log(ViscosityUrbain))/G);
       
    else %Method 1
        
        ElCond_Giordano = exp((-0.08-log(ViscosityGiordano))/1.18);
        ElCond_Optical  = exp((-0.08-log(ViscosityOptical))/1.18);
        ElCond_Riboud = exp((-0.08-log(ViscosityRiboud))/1.18);
        ElCond_Urbain = exp((-0.08-log(ViscosityUrbain))/1.18);
    end

end
