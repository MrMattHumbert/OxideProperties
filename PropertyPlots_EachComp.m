%Properties Plots
load('Prop_Apollo.mat')



%% Make some nice plots!

% Plot for mass of FeO and SiO2 over time
figure

yyaxis left
plot(x,luna.comp_mass(5,:),'-',x,luna.comp_mass(1,:),'--',x,apollo.comp_mass(5,:),':',x,apollo.comp_mass(1,:),'-.');
ylabel('Mass of Species / g/100g');
yyaxis right
plot(x,luna.Lam_Corr,'-',x,apollo.Lam_Corr,'--')
ylabel('Q')
xlabel('Mass of O_{2} Recovered /g');
title('Concentration of FeO and SiO_{2} During Electrolysis at 1800\circC','FontWeight','Normal');
grid on
legend('FeO Luna', 'SiO_{2} Luna','FeO Apollo','SiO_{2} Apollo','Optical Bacisity Luna','Optical Bacisity Apollo','Location','northeast');

% figure
% yyaxis right
%plot(x,luna.Q,'-',x,apollo.Q,'--')
% plot(x,luna.Lam_Corr,'-.',x,apollo.Lam_Corr,':')
% xlabel('Mass of O_{2} Recovered /g');
% ylabel('Optical Basicity')
% legend('Optical Bacisity Luna','Optical Bacisity Apollo')
% 'Q Luna','Q Apollo',





%%
% Melting Point
figure
plot(x,luna.T_liq-273,'-',x,apollo.T_liq-273,'-.'); 
ylabel('Temperature/\circC');
yyaxis right
plot(x,luna.Lam_Corr,'--',x,apollo.Lam_Corr,':')
ylabel('Optical Basicity')
xlabel('Mass of O_{2} Recovered /g');
title('Bath Melting Temperatures and Optical Basicity');
legend('Liquidus Luna','Liquidus Apollo','Basicity Luna','Basicity Apollo','Location','northwest'); %,'Glass Transition');

%%
% Heat Capacity
figure
subplot(1,2,1)
plot(x,luna.Cp_T_Steb,'-',x,luna.Cp_T_L_Steb,'-.',x,luna.Cp_T_Mills,'--'); %,x,Cp_T_NIST,':')
ylabel('Heat Capacity /J-mol^{-1}');
xlabel('Mass of O_{2} Recovered /g');
title('Heat Capacity During Electrolysis at 1800\circC');
grid on
legend('Stebbins','Stebbins Liquid','Mills Liquid','Location','northwest') %,'NIST Liquid'

subplot(1,2,2)
plot(x,apollo.Cp_T_Steb,'-',x,apollo.Cp_T_L_Steb,'-.',x,apollo.Cp_T_Mills,'--'); %,x,Cp_T_NIST,':')
ylabel('Heat Capacity /J-mol^{-1}');
xlabel('Mass of O_{2} Recovered /g');
title('Heat Capacity During Electrolysis at 1800\circC');
grid on
legend('Stebbins','Stebbins Liquid','Mills Liquid','Location','northwest') %,'NIST Liquid'
%%
% Density
figure
plot(x,luna.DensityMills,'-',x,luna.DensityKeen/1000,'--',x,luna.DensityStebbins,x,apollo.DensityMills,':',x,apollo.DensityKeen/1000,'-.',x,apollo.DensityStebbins)
ylabel('Density of Electrolyte/g-cm^{-3}');
xlabel('Mass of O_{2} Recovered /g');
title('Density During Electrolysis at 1800\circC');
grid on
legend('Luna from Mills', 'Luna from Keen','Luna from Stebbins','Apollo from Mills','Apollo from Keen','Apollo from Stebbins','Location','north')
%%
% Resistivity/Conductivity Luna sample
figure
subplot(1,2,2)
plot(x,(1./luna.ElCondOpt),'--',x,(1./luna.ElCondRib),'-.',x,(1./luna.ElCondUrb),':') %,x,(1./ElCondZhang),'-')
ylabel('Resistivity of Electrolyte/\Omega-cm');
xlabel('Mass of O_{2} Recovered /g');
title('Resistivity During Electrolysis at 1800\circC (Luna)');
grid on
legend('Optical','Riboud','Urbain','Location','northeast')

% Resistivity/Conductivity Apollo Sample
%figure
subplot(1,2,1)
plot(x,(1./apollo.ElCondOpt),'--',x,(1./apollo.ElCondRib),'-.',x,(1./apollo.ElCondUrb),':') %,x,(1./ElCondZhang),'-')
ylabel('Resistivity of Electrolyte/\Omega-cm');
xlabel('Mass of O_{2} Recovered /g');
title('Resistivity During Electrolysis at 1800\circC (Apollo)');
grid on
legend('Optical','Riboud','Urbain','Location','northeast')
% Plot Points


%%
% Viscosity Luna 
figure
subplot(1,2,1);
plot(x,luna.ViscosityOptical,'-.',x,luna.ViscosityRiboud,'--',x,luna.ViscosityUrbain,':')
ylabel('Viscosity of Electrolyte / dPas');
xlabel('Mass of O_{2} Recovered / g');
title('Viscosity During Electrolysis at 1800\circC (Luna)');
grid on
legend('Optical','Riboud','Urbain','Location','northeast')


% Viscosity Apollo
%figure
subplot(1,2,2)
plot(x,apollo.ViscosityOptical,'-.',x,apollo.ViscosityRiboud,'--',x,apollo.ViscosityUrbain,':')
ylabel('Viscosity of Electrolyte / dPas');
xlabel('Mass of O_{2} Recovered / g');
title('Viscosity During Electrolysis at 1800\circC (Apollo)');
grid on
legend('Optical','Riboud','Urbain','Location','northeast')
%%
%Viscosity all models - Apollo
figure 
subplot(1,2,1)
plot(x(1:4:end),viscosity_apollo(:,[1,2,4,5,6]))
ylabel('Viscosity of Electrolyte / dPas');
xlabel('Mass of O_{2} Recovered / g');
title('Viscosity During Electrolysis at 1800\circC (Apollo)');
grid on
legend('FactSAGE Glasses','FactSAGE Melts','Optical','Riboud','Urbain','Location','northeast')

%Viscosity all models - Luna
subplot(1,2,2)
plot(x(1:4:end),viscosity_luna(:,[1,2,4,5,6]))
ylabel('Viscosity of Electrolyte / dPas');
xlabel('Mass of O_{2} Recovered / g');
title('Viscosity During Electrolysis at 1800\circC (Luna)');
grid on
legend('FactSAGE Glasses','FactSAGE Melts','Optical','Riboud','Urbain','Location','northeast')

%%
% Surface Tension
figure
%plot(x,SurfTensMills0,'-',x,SurfTensMills1,'--',x,SurfTensMills1b,'--',x,SurfTensMills2,':',x,SurfTensMills2b,'-.')
plot(x,luna.SurfTensMills2b,'-',x,apollo.SurfTensMills2b,'--')
ylabel('Surface Tension / N-m^{-1}');
xlabel('Mass of O_{2} Recovered /g');
title('Surface Tension During Electrolysis at 1800\circC');
grid on
legend('Luna','Apollo','Location','northwest')
%%
% Thermal Conductivity Luna
figure
subplot(1,2,1)
plot(x,luna.ThermalCondMills,'-',x,luna.ThermalCondRiboud_liq,'-.')
ylabel('Thermal Conductivity / W-m^{-1}K^{-1}');
xlabel('Mass of O_{2} Recovered /g');
title('Thermal Conductivity During Electrolysis at 1800\circC (Luna)');
grid on
legend('Mills','Riboud Liquid','Location','northeast')

% Thermal Conductivity Apollo
%figure
subplot(1,2,2)
plot(x,apollo.ThermalCondMills,'-',x,apollo.ThermalCondRiboud_liq,'-.')
ylabel('Thermal Conductivity / W-m^{-1}K^{-1}');
xlabel('Mass of O_{2} Recovered /g');
title('Thermal Conductivity During Electrolysis at 1800\circC (Apollo)');
grid on
legend('Mills','Riboud Liquid','Location','northeast')
