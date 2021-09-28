%{
Ref: Mills, K. C., Yuan, L., Li, Z., & Zhang, G. (2013). Estimating viscosities, electrical & thermal conductivities of slags. High Temperatures - High Pressures, 42(3), 237–256. http://search.ebscohost.com/login.aspx?direct=true&profile=ehost&scope=site&authtype=crawler&jrnl=00181544&AN=89590983&h=qj5veiASLr2Of2c538m29TkZsYPxn3glJbjuvMXfZhPUEOO5reg67DyLchXwFbCLuQPn2NRoZyLl%2FNlos9wCpQ%3D%3D&crl=c&casa_token=jDqO98RZRLQAAAAA:x9AcP

Notes: 

Output in Ohm-cm
%}

function [ResistivityZhang] = Res_Zhang(T,XSiO2,XTiO2,XAl2O3,XCr2O3,XFeO,XMgO,XMnO,XCaO,XNa2O,XK2O,XLi2O,XCaF2,XZrO2,XB2O3,XCrO,XNiO,XFe2O3,XBaO,XSrO)

% Call supporting Functions
% X_mo =  XFeO+...    %FeO
%     XMgO+...    %MgO
%     XMnO+...    %MnO
%     XTiO2*2+...  %TiO2
%     XCaO+...    %CaO
%     XCr2O3*0.45*4;%Cr2O3
% X_m2o =  XNa2O+...   %Na2O
%     XK2O   ;  %K2O
% 
% 
% f_mo = X_mo/(X_mo+X_m2o);
% f_as = XAl2O3/(XAl2O3+XSiO2);
% f_nb = 0.4+0.8*XSiO2;
% 
% Qf = 4-(2*(X_mo+X_m2o+3*f_nb*XFe2O3 - XAl2O3-(1-f_nb)*XFe2O3)/(XSiO2+2*XAl2O3+2*(1-f_nb)*XFe2O3));
[~,Qf,~,~] = Q(T,XSiO2,XTiO2,XAl2O3,XCr2O3,XFeO,XMgO,XMnO,XCaO,XNa2O,XK2O,XLi2O,XCaF2,XZrO2,XB2O3,XCrO,XNiO,XFe2O3,XBaO,XSrO);

% Model

    if XAl2O3>0
        del_n = -XAl2O3;
        
        if Qf<=2.8 %eq 34&35 %Crystalline 
            lnR_1900 = (-5.37 + 2.674*exp(Qf/3.43)+ 3.5e-7*exp(Qf/0.239))+(3.104+5.75*del_n-25.8*del_n^2 + 31*del_n^3);
            B_R = (7.34+0.553*exp(Qf/1.0123))+(30.03-60.5*del_n+29.34*del_n^2);
            ResistivityZhang = exp(lnR_1900 + (1000*B_R/T)-(1000*B_R/1900));
        else % eq 37&38 %Glassy
            lnR_1900 = (-0.55+0.635*Qf) + (3.104+5.75*del_n-25.8*del_n^2 + 31*del_n^3);
            B_R  = (13.0+2.7*Qf) + (30.03-60.5*del_n+29.34*del_n^2);
            ResistivityZhang = exp(lnR_1900 + (1000*B_R/T)-(1000*B_R/1900));
        end
    else
        del_n = X_m2o;
        if del_n>=0.1  
            lnR_1900 = (-4.37 + 2.674*exp(Qf/3.43)+ 3.5e-7*exp(Qf/0.239))+(3.104+5.75*del_n-25.8*del_n^2 + 31*del_n^3);
            B_R = (7.34+0.553*exp(Qf/1.0123))+(30.03-60.5*del_n+29.34*del_n^2);
            ResistivityZhang = exp(lnR_1900 + (1000*B_R/T)-(1000*B_R/1900));
        else 
            lnR_1900 = (-4.37 + 2.674*exp(Qf/3.43) + 3.5e-7*exp(Qf/0.239))+33.6*del_n;
            B_R = (7.34+0.553*exp(Qf/1.0123))+ 240*del_n;
            ResistivityZhang = exp(lnR_1900 + (1000*B_R/T)-(1000*B_R/1900));
        end
    end

end
