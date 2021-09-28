%{
Ref: Mills, K. C., Yuan, L., Li, Z., & Zhang, G. (2013). Estimating viscosities, electrical & thermal conductivities of slags. High Temperatures - High Pressures, 42(3), 237–256. http://search.ebscohost.com/login.aspx?direct=true&profile=ehost&scope=site&authtype=crawler&jrnl=00181544&AN=89590983&h=qj5veiASLr2Of2c538m29TkZsYPxn3glJbjuvMXfZhPUEOO5reg67DyLchXwFbCLuQPn2NRoZyLl%2FNlos9wCpQ%3D%3D&crl=c&casa_token=jDqO98RZRLQAAAAA:x9AcP

Notes: 

Output in W/(m-K)
%}

function [ThermCond_Zhang] =  ThermCond_Zhang(T,XSiO2,XTiO2,XAl2O3,XCr2O3,XFeO,XMgO,XMnO,XCaO,XNa2O,XK2O,XLi2O,XCaF2,XZrO2,XB2O3,XCrO,XNiO,XFe2O3,XBaO,XSrO)

%Supporting Functions
[~,T_Crit,T_Liq,~] = Temps(T,XSiO2,XTiO2,XAl2O3,XCr2O3,XFeO,XMgO,XMnO,XCaO,XNa2O,XK2O,XLi2O,XCaF2,XZrO2,XB2O3,XCrO,XNiO,XFe2O3,XBaO,XSrO);

feed_mol_frac = [XSiO2,XTiO2,XAl2O3,XCr2O3,XFeO,XMgO,XMnO,XCaO,XNa2O,XK2O,XLi2O,XCaF2,XZrO2,XB2O3,XCrO,XNiO,XFe2O3,XBaO,XSrO]';

X_mo =  XFeO+...    %FeO
    XMgO+...    %MgO
    XMnO+...    %MnO
    XTiO2*2+...  %TiO2
    XCaO+...    %CaO
    XCr2O3*0.45*4;%Cr2O3
X_m2o =  XNa2O+...   %Na2O
    XK2O   ;  %K2O

f_mo = X_mo/(X_mo+X_m2o);
f_as = XAl2O3/(XAl2O3+XSiO2);
f_nb = 0.4+0.8*XSiO2;

Qf = 4-(2*(X_mo+X_m2o+3*f_nb*XFe2O3 - XAl2O3-(1-f_nb)*XFe2O3)/(XSiO2+2*XAl2O3+2*(1-f_nb)*XFe2O3));


%         SiO2, TiO2,Al2O3,Cr2O3,  FeO,   MgO, MnO, CaO,Na2O, K2O,Li2O,CaF2, ZrO2,B2O3, CrO, NiO,Fe2O3,BaO,  SrO
r       =[0.40,0.605,0.535,0.615,0.745, 0.72, 0.67, 1.0,1.02,1.38,0.76,   1, 0.72,0.27,0.73,0.69,0.55,1.35, 1.18]'; %10^-10 m
z       =[   2,    4,   3,     3,    2,    2,    2,   2,   1,   1,   1,   2,    4,   3,   2,   2,   3,   2,   2 ]';
z_r2    = z./(r.^2);
rm_rCa3 = r.^3; 
      Zmean = [0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1]'; % rectify with table in Zhang
rm_rca_mean = sum(feed_mol_frac.*rm_rCa3.*Zmean)/sum(feed_mol_frac.*Zmean);
z_r2_mean = sum(z_r2.*Zmean.*feed_mol_frac);

% Model
    k_298 = 7.94 + 7.219*exp(Qf/199.3)+0.109*exp(Qf/3.634)+0.612*(X_mo+X_m2o-XAl2O3)*z_r2_mean;
    k_Tcrit = -0.851+0.987*exp(Qf/1.006)-0.80*exp(Qf/0.9646)+1.621*(X_mo+X_m2o-XAl2O3)*z_r2_mean;
    k_Tliq = 0.139+3.65e-5*exp(Qf/0.3421);
    k_1773 = 0.056+8.149e-10*exp(Qf/0.179);
    
ThermCond_Zhang = polyval(polyfit([298, T_Crit, T_Liq, 1773],[k_298,k_Tcrit,k_Tliq,k_1773],2),T);

end
