%{
Ref: Mills, K. C. (2011). The Estimation Of Slag Properties. Southern African Pyrometallurgy 2011 International Conference, March, 1–52. http://www.pyrometallurgy.co.za/KenMills/KenMills.pdf 
Ref: Mills, K. C., Karagadde, S., Lee, P. D., Yuan, L., & Shahbazian, F. (2016). Calculation of physical properties for use in models of continuous casting process-Part 1: Mould slags. ISIJ International, 56(2). https://doi.org/10.2355/isijinternational.ISIJINT-2015-364

Notes: Reproduced from Mills

Output in mN/m
%}

function [SurfTensMills2, SurfTensMills2b] = SurfTens_Method2(T,XSiO2,XTiO2,XAl2O3,XCr2O3,XFeO,XMgO,XMnO,XCaO,XNa2O,XK2O,XLi2O,XCaF2,XZrO2,XB2O3,XCrO,XNiO,XFe2O3,XBaO,XSrO)
    %           SiO2, TiO2, Al2O3, Cr2O3, FeO,  MgO, MnO,   CaO, Na2O,  K2O, Li2O, CaF2, ZrO2,B2O3,CrO,NiO,Fe2O3, BaO, SrO
    
    feed_mol_frac = [XSiO2,XTiO2,XAl2O3,XCr2O3,  XFeO,  XMgO,XMnO,  XCaO,XNa2O, XK2O,XLi2O,XCaF2, XZrO2,XB2O3,XCrO,XNiO, XFe2O3,  XBaO,  XSrO]';
    gam  =[      260,  360,   655,   800, 645,  635, 645,   625,  295,  160,  300,  290,  400, 110,645,645,  300, 560, 600]';
    dGam_dT =[ 0.031,-0.15,-0.177, 0.125, 0.1,-0.13,-0.1,-0.094,-0.11,-0.11,-0.11,-0.07,-0.15,0.04,  0,  0,-0.05,-0.1,-0.1]';
    
    % Method 2 
    feed_mol_frac_method2 = feed_mol_frac;
    if (0.12-XB2O3)<0
        feed_mol_frac_method2(14)=abs(0.12-XB2O3);
        X_gam = feed_mol_frac_method2.*gam;
        X_gam(14)= X_gam(14)+265.3*(0.12) -155.3;
    elseif (0.12-XB2O3-XK2O)<0
        feed_mol_frac_method2(10)=abs(0.12-XB2O3-XK2O);
        X_gam = feed_mol_frac_method2.*gam;
        if 0<(0.12-XB2O3) && (0.12-XB2O3)<0.115 
            X_gam(10) =X_gam(10)+0.8-1388*(0.12-XB2O3)+6723*(0.12-XB2O3)^2;
        else
            X_gam(10)=X_gam(10)+254.5*(0.12-XB2O3)-94.5;
        end    
    elseif  (0.12-XB2O3-XK2O-XNa2O)<0
        feed_mol_frac_method2(9)=abs(0.12-XB2O3-XK2O-XNa2O);
        X_gam = feed_mol_frac_method2.*gam;
        if 0<(0.12-XB2O3-XK2O) && (0.12-XB2O3-XK2O)<0.115   
            X_gam(9) =X_gam(9)+ 0.8-1388*(0.12-XB2O3-XK2O)+6723*(0.12-XB2O3-XK2O)^2;
        else
            X_gam(9)=X_gam(9)-412.9*(0.12-XB2O3-XK2O)-115.9;
        end   
    elseif (0.12-XB2O3-XK2O-XNa2O-XCaF2)<0
        feed_mol_frac_method2(12)=abs(0.12-XB2O3-XK2O-XNa2O-XCaF2);
        X_gam = feed_mol_frac_method2.*gam;
        if 0<(0.12-XB2O3-XK2O-XNa2O) && (0.12-XB2O3-XK2O-XNa2O)<0.13 
            X_gam(12) =X_gam(12)-2-934*(0.12-XB2O3-XK2O-XNa2O)+4769*(0.12-XB2O3-XK2O-XNa2O)^2;
        else
            X_gam(12)=X_gam(12)+382.5*(0.12-XB2O3-XK2O-XNa2O)-92.5;
        end 
    elseif (0.12-XB2O3-XK2O-XNa2O-XCaF2-XFe2O3)<0
        feed_mol_frac_method2(12)=abs(0.12-XB2O3-XK2O-XNa2O-XCaF2-XFe2O3);
        X_gam = feed_mol_frac_method2.*gam;
        if 0<(0.12-XB2O3-XK2O-XNa2O-XCaF2) && (0.12-XB2O3-XK2O-XNa2O-XCaF2)<0.125 
            X_gam(17) =X_gam(17) -2972*(0.12-XB2O3-XK2O-XNa2O-XCaF2)+14312*(0.12-XB2O3-XK2O-XNa2O-XCaF2)^2-3.7;
        else
            X_gam(17)=X_gam(17)+516.2*(0.12-XB2O3-XK2O-XNa2O-XCaF2)-216.2;
        end
    elseif (0.12-XB2O3-XK2O-XNa2O-XCaF2-XFe2O3-XCr2O3)<0
        feed_mol_frac_method2(4)=abs(0.12-XB2O3-XK2O-XNa2O-XCaF2-XFe2O3-XCr2O3);
        X_gam = feed_mol_frac_method2.*gam;
        if 0<(0.12-XB2O3-XK2O-XNa2O-XCaF2-XFe2O3) && (0.12-XB2O3-XK2O-XNa2O-XCaF2-XFe2O3)<0.05 
            X_gam(4)=X_gam(4) -1248*(0.12-XB2O3-XK2O-XNa2O-XCaF2-XFe2O3)+8735*(0.12-XB2O3-XK2O-XNa2O-XCaF2-XFe2O3)^2;
        else
            X_gam(4)=X_gam(4)+884.2*(0.12-XB2O3-XK2O-XNa2O-XCaF2-XFe2O3)-84.2;
        end
    else %Method 1 and 2 are the same
        X_gam = feed_mol_frac.*gam;
        if XNa2O<eps % for Na2O
            X_gam(9) = 0;
        elseif 0<XNa2O && XNa2O<0.115
            X_gam(9) = 0.8-1388*XNa2O+6723*(XNa2O)^2;
        else
            X_gam(9)=-412.9*XNa2O-115.9;
        end
        
        if XK2O<eps % for K2O
            X_gam(10) = 0;
        elseif 0<XK2O && XK2O<0.115
            X_gam(10) = 0.8-1388*XK2O+6723*(XK2O)^2;
        else
            X_gam(10)=254.5*XK2O-94.5;
        end
        
        if XCaF2<eps % for CaF2
            X_gam(12) = 0;
        elseif 0<XCaF2 && XCaF2<0.13
            X_gam(12) = -2-934*XCaF2+4769*(XCaF2)^2;
        else
            X_gam(12)=382.5*XCaF2-92.5;
        end
        
        if XB2O3<eps % for B2O3
            X_gam(14) = 0;
        elseif 0<XB2O3 && XB2O3<0.1
            X_gam(14) = -5.2-3454*XB2O3+22178*(XB2O3)^2;
        else
            X_gam(14)= 265.3*XB2O3 -155.3;
        end
        
        if XCr2O3<eps % for Cr2O3
            X_gam(4) = 0;
        elseif 0<XCr2O3 && XCr2O3<0.05
            X_gam(4) = -1248*XCr2O3+8735*XCr2O3^2;
        else
            X_gam(4)=884.2*XCr2O3-84.2;
        end
        
        if XFe2O3<eps % for Fe2O3
            X_gam(17) = 0;
        elseif 0<XFe2O3 && XFe2O3<0.125
            X_gam(17) = -2972*XFe2O3+14312*XFe2O3^2-3.7;
        else
            X_gam(17)=516.2*XFe2O3-216.2;
        end
    end
    SurfTensMills2 = sum(X_gam)-0.15*(T-1773); % T in K
    SurfTensMills2b = sum(X_gam)-sum(feed_mol_frac.*dGam_dT)*(T-1773);
end