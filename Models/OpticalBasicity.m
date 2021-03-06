%{
Ref: Mills, K. C. (2011). The Estimation Of Slag Properties. Southern African Pyrometallurgy 2011 International Conference, March, 1–52. http://www.pyrometallurgy.co.za/KenMills/KenMills.pdf
Ref: Mills, K. C., Karagadde, S., Lee, P. D., Yuan, L., & Shahbazian, F. (2016). Calculation of physical properties for use in models of continuous casting process-Part 1: Mould slags. ISIJ International, 56(2). https://doi.org/10.2355/isijinternational.ISIJINT-2015-364

Ref: 
Ref: 


Notes: 

Output is dimensionless
%}

function [Lam_Corr] =  OpticalBasicity(T,XSiO2,XTiO2,XAl2O3,XCr2O3,XFeO,XMgO,XMnO,XCaO,XNa2O,XK2O,XLi2O,XCaF2,XZrO2,XB2O3,XCrO,XNiO,XFe2O3,XBaO,XSrO)

XNB = XCaO+XMgO+XNa2O+XK2O+XLi2O+XFeO+XMnO+XCaF2+XCr2O3+XCrO+XNiO+XFe2O3+XBaO;
f = (XNB-XAl2O3)/XNB;
modX_n = [  XSiO2*2,... 
            XTiO2*2,...
            XAl2O3*3,... 
            XCr2O3*f*3,... 
            XFeO*f,... 
            XMgO*f,... 
            XMnO*f,...
            XCaO*f,...
            XNa2O*f,...
            XK2O*f,... 
            XLi2O*f,...
            XCaF2,...
            XZrO2*2,...
            XB2O3*3,...
            XCrO*f,... 
            XNiO*f,... 
            XFe2O3*f*3,...
            XBaO*f,... 
            XSrO*f];
modX_nLam = [  XSiO2*2*0.48,... 
            XTiO2*2*0.61,...
            XAl2O3*3*0.6,... 
            XCr2O3*f*3*0.75,... 
            XFeO*f,... 
            XMgO*f*0.78,... 
            XMnO*f,...
            XCaO*f,...
            XNa2O*f*1.15,...
            XK2O*f*1.4,... 
            XLi2O*f,...
            XCaF2*1.2,...
            XZrO2*2*0.61,...
            XB2O3*3*0.42,...
            XCrO*f,... 
            XNiO*f,... 
            XFe2O3*f*3*0.75,...
            XBaO*f*1.15,... 
            XSrO*f*1.1];
Lam_Corr = sum(modX_nLam)/sum(modX_n);

end