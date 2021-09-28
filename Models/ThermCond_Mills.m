%{
Ref: Mills, K. C. (2011). The Estimation Of Slag Properties. Southern African Pyrometallurgy 2011 International Conference, March, 1–52. http://www.pyrometallurgy.co.za/KenMills/KenMills.pdf
Ref: Mills KC, Karagadde S, Lee PD, Yuan L, Shahbazian F. Calculation of physical properties for use in models of continuous casting process-part 1: mould slags. ISIJ international. 2016 Feb 15;56(2):264-73.

Notes: Q model should not be used for Slags with Q>3.3 or Q<2
        Karagadde provides a more complete model 

Output in W/(m-K)
%}

function [ThermCond_Riboud,ThermCond_Mills,ThermCond_Karagadde,ThermCond_Urbain] = ThermCond_Mills(T,XSiO2,XTiO2,XAl2O3,XCr2O3,XFeO,XMgO,XMnO,XCaO,XNa2O,XK2O,XLi2O,XCaF2,XZrO2,XB2O3,XCrO,XNiO,XFe2O3,XBaO,XSrO)

%Supporting Functions
[~,Qf,~,~] = Q(T,XSiO2,XTiO2,XAl2O3,XCr2O3,XFeO,XMgO,XMnO,XCaO,XNa2O,XK2O,XLi2O,XCaF2,XZrO2,XB2O3,XCrO,XNiO,XFe2O3,XBaO,XSrO);
ViscosityRiboud = Vis_Riboud(T,XSiO2,XTiO2,XAl2O3,XCr2O3,XFeO,XMgO,XMnO,XCaO,XNa2O,XK2O,XLi2O,XCaF2,XZrO2,XB2O3,XCrO,XNiO,XFe2O3,XBaO,XSrO);
ViscosityUrbain = Vis_Urbain(T,XSiO2,XTiO2,XAl2O3,XCr2O3,XFeO,XMgO,XMnO,XCaO,XNa2O,XK2O,XLi2O,XCaF2,XZrO2,XB2O3,XCrO,XNiO,XFe2O3,XBaO,XSrO);
[~,T_Critical,T_Liq,~] = Temps(T,XSiO2,XTiO2,XAl2O3,XCr2O3,XFeO,XMgO,XMnO,XCaO,XNa2O,XK2O,XLi2O,XCaF2,XZrO2,XB2O3,XCrO,XNiO,XFe2O3,XBaO,XSrO);
DensityMills = Rho_Mills(T,XSiO2,XTiO2,XAl2O3,XCr2O3,XFeO,XMgO,XMnO,XCaO,XNa2O,XK2O,XLi2O,XCaF2,XZrO2,XB2O3,XCrO,XNiO,XFe2O3,XBaO,XSrO);
HeatCapacity = Cp_NIST(T,XSiO2,XTiO2,XAl2O3,XCr2O3,XFeO,XMgO,XMnO,XCaO,XNa2O,XK2O,XLi2O,XCaF2,XZrO2,XB2O3,XCrO,XNiO,XFe2O3,XBaO,XSrO);
feed_mol_frac =  [XSiO2, XTiO2, XAl2O3, XCr2O3,   XFeO,   XMgO, XMnO,   XCaO, XNa2O,  XK2O, XLi2O, XCaF2,  XZrO2, XB2O3, XCrO,  XNiO, XFe2O3,   XBaO,   XSrO]';
Mol_Weight = [  60, 79.9, 101.9,   152, 71.85,  40.3,  71, 56.08,   62, 94.2, 29.8,   78, 123.2, 69.8,  68, 74.7, 159.7, 153.3, 103.6]';


% Riboud Model
ThermCond_Riboud = exp(-2.178+(0.2821*log(ViscosityRiboud)));
% Urbain Model
ThermCond_Urbain = exp(-2.178+(0.2821*log(ViscosityUrbain)));

% Karagadde for molten slags
ThermCond_Karagadde = 0.139+3.65e-5*exp(Qf/0.3421);

if T<T_Liq
    % Solid Slag Model
    if Qf<=2.8 %Crystalline
        if T<T_Critical
            k_298 = exp(-0.424+0.00002*exp(Qf/0.299)+3.2*XLi2O);
            k_crit = exp(-0.435+0.00005*exp(Qf/0.332)+3*XLi2O);
            ThermCond_Mills= k_298-(T-298)*(k_crit-k_298)/(T_Critical-298);
        else % T_Critical < T < T_Liq
            k_crit = exp(-0.435+0.00005*exp(Qf/0.332)+3*XLi2O);
            k_liq = exp(-1.914+0.00037*exp(Qf/0.402));
            ThermCond_Mills = k_crit-(T-T_Critical)*(k_liq-k_crit)/(T-T_Critical);
        end
    else %Glassy
        ThermCond_Mills = 4e-1 * DensityMills/sum(feed_mol_frac.*Mol_Weight) * HeatCapacity;
    end
else % Molten Slag
        % Mills Q model
        dk_dT =(5e-4)*exp(Qf*0.5551);
        ThermCond_Mills = exp(-1.914+0.00037*exp(Qf/0.402))+dk_dT*(T-T_Liq);
end

    
end

