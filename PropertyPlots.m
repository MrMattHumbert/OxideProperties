%Properties Plots
clc;
clear all;

load('Prop_Apollo.mat')
load('Prop_Luna.mat')


x_A = Apollo.Trajectory.MaxO2_g*linspace(0,1,size(Apollo.Trajectory.MolEvolution,2)); 
x_L = Luna.Trajectory.MaxO2_g*linspace(0,1,size(Luna.Trajectory.MolEvolution,2));
% From Main script T = [1800,1850,1900,1950]+273 ; % in Kelvin can be single point or vector
T = 2; % Which Temperature point should be used for plotting
%% Make some nice plots!

% Plot for mass of FeO and SiO2 over time MassEvolution is a 19 component
% vector of all the species
figure

yyaxis left
plot(x_A,Apollo.Trajectory.MassEvolution(5,:),'-',x_A,Apollo.Trajectory.MassEvolution(1,:),'--',x_L,Luna.Trajectory.MassEvolution(5,:),'-.',x_L,Luna.Trajectory.MassEvolution(1,:),':');
ylabel('Mass of Species / g/100g');
yyaxis right
plot(x_A,Apollo.OpticalBasicity(T,:),x_L,Luna.OpticalBasicity(T,:)) % Basicity doesn't change with temp but good to keep consistant
ylabel('Corrected Optical Basicity')
xlabel('Mass of O_{2} Recovered / g/100g');
title('Concentration of FeO and SiO_{2} During Electrolysis','FontWeight','Normal');
grid on
legend('FeO Apollo', 'SiO_{2} Apollo','FeO Luna','SiO_{2} Luna','Basicity Apollo','Basicity Luna','Location','North');

saveas(gcf,'Plots/ConcentrationBasicity.svg')


%%
% Melting Point
figure
plot(x_A,Apollo.Temperatures.T_liq(T,:)-273,'-',x_L,Luna.Temperatures.T_liq(T,:)-273,'-.'); 
% hold on;
% plot(x_A,Apollo.Temperatures.T_crit(1,:)-273,'-',x_L,Luna.Temperatures.T_crit(1,:)-273,'-.'); 
% plot(x_A,Apollo.Temperatures.T_g(1,:)-273,'-',x_L,Luna.Temperatures.T_g(1,:)-273,'-.'); 
ylabel('Temperature/\circC');
yyaxis right
plot(x_A,Apollo.Q(T,:),'--',x_L,Luna.Q(T,:),':') % Q doesn't change with temp, but good to keep consistant
ylabel('Q')
xlabel('Mass of O_{2} Recovered / g/100g');
title('Bath Liquidus Temperatures and Degree of Polymerization','FontWeight','Normal');
legend('Liquidus Apollo','Liquidus Luna','Q Apollo','Q Luna','Location','NorthEast'); 

saveas(gcf, 'Plots/MeltingPoint.svg')
%%
% Heat Capacity
figure

subplot(1,2,1)
plot(x_L,Luna.HeatCapacity.Stebbins(T,:),'-',x_L,Luna.HeatCapacity.StebbinsLiq(T,:),'-.',x_L,Luna.HeatCapacity.Mills(T,:),'--',x_L,Luna.HeatCapacity.NIST(T,:),':');
ylabel('Heat Capacity /J-mol^{-1}');
xlabel('Mass of O_{2} Recovered / g/100g');
title('Heat Capacity at 1850\circC (Luna)','FontWeight','Normal');
grid on
legend('Stebbins','Stebbins Liquid','Mills','NIST','Location','northwest') 

subplot(1,2,2)
plot(x_A,Apollo.HeatCapacity.Stebbins(T,:),'-',x_A,Apollo.HeatCapacity.StebbinsLiq(T,:),'-.',x_A,Apollo.HeatCapacity.Mills(T,:),'--',x_A,Apollo.HeatCapacity.NIST(T,:),':');
ylabel('Heat Capacity /J-mol^{-1}');
xlabel('Mass of O_{2} Recovered / g/100g');
title('Heat Capacity at 1850\circC (Apollo)','FontWeight','Normal');
grid on
legend('Stebbins','Stebbins Liquid','Mills','NIST','Location','northwest') 

saveas(gcf,'Plots/HeatCapacity.svg')

%%
% Density
figure
plot(x_L,Luna.Density.Mills(T,:),'-',x_L,Luna.Density.Keen(T,:),'--',x_L,Luna.Density.Stebbins(T,:),'-.',x_A,Apollo.Density.Mills(T,:),':',x_A,Apollo.Density.Keen(T,:),'-.',x_A,Apollo.Density.Stebbins(T,:))
ylabel('Density of Electrolyte/g-cm^{-3}');
xlabel('Mass of O_{2} Recovered / g/100g');
title('Density at 1850\circC','FontWeight','Normal');
grid on
legend('Luna from Mills', 'Luna from Keen','Luna from Stebbins','Apollo from Mills','Apollo from Keen','Apollo from Stebbins','Location','north')

saveas(gcf,'Plots/Density.svg')
%%
% Resistivity/Conductivity Luna sample
figure

subplot(1,2,1)
plot(x_L,(Luna.Flow.Res_Zhang(T,:)),'-', x_L,(1./Luna.Flow.ElCond_Giordano(T,:)),'--',x_L,(1./Luna.Flow.ElCond_Optical(T,:)),':',x_L,(1./Luna.Flow.ElCond_Riboud(T,:)),':',x_L,(1./Luna.Flow.ElCond_Urbain(T,:)),':'); 
ylabel('Resistivity of Electrolyte/\Omega-cm');
xlabel('Mass of O_{2} Recovered / g/100g');
title('Resistivity at 1850\circC (Luna)','FontWeight','Normal');
grid on
legend('Zhang','Giordano','Optical','Riboud','Urbain','Location','northeast')

% Resistivity/Conductivity Apollo Sample
%figure
subplot(1,2,2)
plot(x_A,(Apollo.Flow.Res_Zhang(T,:)),'-',x_A,(1./Apollo.Flow.ElCond_Giordano(T,:)),'--',x_A,(1./Apollo.Flow.ElCond_Optical(T,:)),':',x_A,(1./Apollo.Flow.ElCond_Riboud(T,:)),':',x_A,(1./Apollo.Flow.ElCond_Urbain(T,:)),':');
ylabel('Resistivity of Electrolyte/\Omega-cm');
xlabel('Mass of O_{2} Recovered / g/100g');
title('Resistivity at 1850\circC (Apollo)','FontWeight','Normal');
grid on
legend('Zhang','Giordano','Optical','Riboud','Urbain','Location','northeast')

saveas(gcf,'Plots/Resistivity.svg')


%%
% Viscosity Luna 
figure

subplot(1,2,1);
plot(x_L,(Luna.Flow.ViscosityZhang(T,:)),'-', x_L,(Luna.Flow.ViscosityGiordanoNBO(T,:)),'--',x_L,(Luna.Flow.ViscosityGiordanoSM(T,:)),'--',x_L,(Luna.Flow.ViscosityOptical(T,:)),':',x_L,(Luna.Flow.ViscosityRiboud(T,:)),':',x_L,(Luna.Flow.ViscosityUrbain(T,:)),':'); 
ylabel('Viscosity of Electrolyte / dPas');
xlabel('Mass of O_{2} Recovered / g/100g');
title('Viscosity at 1850\circC (Luna)','FontWeight','Normal');
grid on
legend('Zhang','Giordano NBO','Giordano SM','Optical','Riboud','Urbain','Location','northeast')


% Viscosity Apollo
%figure
subplot(1,2,2)
plot(x_A,(Apollo.Flow.ViscosityZhang(T,:)),'-', x_A,(Apollo.Flow.ViscosityGiordanoNBO(T,:)),'--',x_A,(Apollo.Flow.ViscosityGiordanoSM(T,:)),'--',x_A,(Apollo.Flow.ViscosityOptical(T,:)),':',x_A,(Apollo.Flow.ViscosityRiboud(T,:)),':',x_A,(Apollo.Flow.ViscosityUrbain(T,:)),':'); 
ylabel('Viscosity of Electrolyte / dPas');
xlabel('Mass of O_{2} Recovered / g/100g');
title('Viscosity at 1850\circC (Apollo)','FontWeight','Normal');
grid on
legend('Zhang','Giordano NBO','Giordano SM','Optical','Riboud','Urbain','Location','northeast')

saveas(gcf,'Plots/Viscosity.svg')

%%
% Surface Tension
figure
%plot(x,SurfTensMills0,'-',x,SurfTensMills1,'--',x,SurfTensMills1b,'--',x,SurfTensMills2,':',x,SurfTensMills2b,'-.')
plot(x_L,Luna.Surface.SurfTens_Method2b(T,:),'-',x_A,Apollo.Surface.SurfTens_Method2b(T,:),'--')
ylabel('Surface Tension / N-m^{-1}');
xlabel('Mass of O_{2} Recovered / g/100g');
title('Surface Tension at 1850\circC','FontWeight','Normal');
grid on
legend('Luna','Apollo','Location','northwest')

saveas(gcf,'Plots/SurfaceTension.svg')

%%
% Thermal Conductivity Luna
figure
subplot(1,2,1)
plot(x_L,Luna.HeatTransfer.ThermalCond_Mills(T,:),'-',x_L,Luna.HeatTransfer.ThermalCond_Riboud(T,:),'-.',x_L,Luna.HeatTransfer.ThermalCond_Zhang(T,:),':',x_L,Luna.HeatTransfer.ThermalCond_Urbain(T,:),'--')
ylabel('Thermal Conductivity / W-m^{-1}K^{-1}');
xlabel('Mass of O_{2} Recovered / g/100g');
title('Thermal Conductivity at 1850\circC (Luna)','FontWeight','Normal');
grid on
legend('Mills','Riboud','Zhang','Urbain','Location','northeast')

% Thermal Conductivity Apollo
%figure
subplot(1,2,2)
plot(x_A,Apollo.HeatTransfer.ThermalCond_Mills(T,:),'-',x_A,Apollo.HeatTransfer.ThermalCond_Riboud(T,:),'-.',x_A,Apollo.HeatTransfer.ThermalCond_Zhang(T,:),':',x_A,Apollo.HeatTransfer.ThermalCond_Urbain(T,:),'--')
ylabel('Thermal Conductivity / W-m^{-1}K^{-1}');
xlabel('Mass of O_{2} Recovered / g/100g');
title('Thermal Conductivity at 1850\circC (Apollo)','FontWeight','Normal');
grid on
legend('Mills','Riboud','Zhang','Urbain','Location','northeast')

saveas(gcf,'Plots/ThermalConductivity.svg')