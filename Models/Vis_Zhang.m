%{
Ref: Mills, K. C., Yuan, L., Li, Z., & Zhang, G. (2013). Estimating viscosities, electrical & thermal conductivities of slags. High Temperatures - High Pressures, 42(3), 237–256. http://search.ebscohost.com/login.aspx?direct=true&profile=ehost&scope=site&authtype=crawler&jrnl=00181544&AN=89590983&h=qj5veiASLr2Of2c538m29TkZsYPxn3glJbjuvMXfZhPUEOO5reg67DyLchXwFbCLuQPn2NRoZyLl%2FNlos9wCpQ%3D%3D&crl=c&casa_token=jDqO98RZRLQAAAAA:x9AcP

Notes: This viscosity model gives very different results to other models.
It was not reviewed by Karagadde

Output in dPa-s
%}

function [ViscosityZhang] =  Vis_Zhang(T,XSiO2,XTiO2,XAl2O3,XCr2O3,XFeO,XMgO,XMnO,XCaO,XNa2O,XK2O,XLi2O,XCaF2,XZrO2,XB2O3,XCrO,XNiO,XFe2O3,XBaO,XSrO)

%Supporting Functions
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
rm_rCa3 = r.^3; 
      Zmean = [0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1]'; % rectify with table in Zhang
rm_rca_mean = (sum(feed_mol_frac.*rm_rCa3.*Zmean)/sum(feed_mol_frac.*Zmean))^3;


%Model
    if Qf<= 2.8
        ln_eta_1900 = -2.161+0.8279*exp(Qf/1.794)+2.19e-16*exp(Qf/0.1036)+0.94*rm_rca_mean;
        B_a_lowQf = f_mo*(-795830.8+8.476e-6*exp(Qf/0.2611)+795841.4*exp(Qf/211587.4))+...
            (1-f_mo)*(-24.6+6.122e-13*exp(Qf/0.1258)+28.465*exp(Qf/7.727))+...
            (-4.76+4.57*rm_rca_mean);
        ln_eta_T = ln_eta_1900+(1000*B_a_lowQf/T)-(1000*B_a_lowQf/1900);
        ViscosityZhang = exp(ln_eta_T);
    else
        ln_eta_1900 = 18.27 - 114.3*f_as + 282.4*f_as^2 - 232*f_as^3; %ea 19
        ln_eta_1900_s = -2.161+0.8279*exp(Qf/1.794)+2.19e-16*exp(Qf/0.1036)+0.94*rm_rca_mean; %eq 11
        f_curve = (ln_eta_1900-1.782)/(18.37-1.782);
        ln_eta_1900as = 1.782+f_curve*(ln_eta_1900_s-1.782)+0.94*rm_rca_mean;
        B_a_s = 62.08-235.2*f_as+544.5*f_as^2-440.2*f_as^3;
        f_curve_B = (B_a_s-16.3)/(62.6-16.3);
        B_star = 16.3 + (f_curve_B*B_a_s-16.3)-4.76+4.57*rm_rca_mean;
        ln_eta_T = ln_eta_1900as + (1000*B_star/T)-(1000*B_star/1900);
        ViscosityZhang = exp(ln_eta_T);
    end



end