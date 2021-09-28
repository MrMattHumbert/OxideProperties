 %%%%%%%%%%%%%%%%%%%%%%%
% Property evolution when removing regolith
%%%%%%%%%%%%%%%%%%%%%%%

function Comp = PropertiesDuringO2Removal(T,Comp)

%Loop over Temperature and Composition

for t = 1:length(T)

    for x = 1:size(Comp.Trajectory.MolFracEvolution,2)
        X = num2cell(Comp.Trajectory.MolFracEvolution(:,x));
        %OpticalBasicity
        Comp.OpticalBasicity(t,x) = OpticalBasicity(T(t),X{:});
        %Q
        [Comp.Q(t,x),Comp.Qf(t,x),Comp.NBO_T(t,x),Comp.NBO_T2(t,x)] = Q(T(t),X{:});
        %T_break,T_crit,T_liq,T_g]
        [Comp.Temperatures.T_break(t,x),Comp.Temperatures.T_crit(t,x),Comp.Temperatures.T_liq(t,x),Comp.Temperatures.T_g(t,x)] = Temps(T(t),X{:});
 
        %Heat Capacities
        Comp.HeatCapacity.Mills(t,x) = Cp_Mills(T(t),X{:});
        [Comp.HeatCapacity.NIST(t,x),~] = Cp_NIST(T(t),X{:});
        [Comp.HeatCapacity.Stebbins(t,x),~,Comp.HeatCapacity.StebbinsLiq(t,x)] = Cp_Stebbins(T(t),X{:});
        
        %Surface
        [Comp.Surface.SurfTens_Method0(t,x),Comp.Surface.SurfTens_Method1(t,x),Comp.Surface.SurfTens_Method1b(t,x)] = SurfTens_Method1(T(t),X{:});
        [Comp.Surface.SurfTens_Method2(t,x),Comp.Surface.SurfTens_Method2b(t,x)] = SurfTens_Method2(T(t),X{:});
        
        %Heat Transfer
        Comp.HeatTransfer.ThermalCond_Zhang(t,x) = ThermCond_Zhang(T(t),X{:});
        [Comp.HeatTransfer.ThermalCond_Riboud(t,x),Comp.HeatTransfer.ThermalCond_Mills(t,x),Comp.HeatTransfer.ThermalCond_Karagadde(t,x),Comp.HeatTransfer.ThermalCond_Urbain(t,x)] = ThermCond_Mills(T(t),X{:});
        
        
        % Viscosity
        [Comp.Flow.ViscosityGiordanoNBO(t,x),Comp.Flow.ViscosityGiordanoSM(t,x)] = Vis_Giordano(T(t),X{:});
        Comp.Flow.ViscosityZhang(t,x) = Vis_Zhang(T(t),X{:});
        Comp.Flow.ViscosityUrbain(t,x) = Vis_Urbain(T(t),X{:});
        Comp.Flow.ViscosityRiboud(t,x) = Vis_Riboud(T(t),X{:});
        Comp.Flow.ViscosityOptical(t,x) = Vis_Optical(T(t),X{:});
        
        %Electrical Conductivity
        Comp.Flow.Res_Zhang(t,x) = Res_Zhang(T(t),X{:});
        [Comp.Flow.ElCond_Giordano(t,x),Comp.Flow.ElCond_Optical(t,x),Comp.Flow.ElCond_Riboud(t,x),Comp.Flow.ElCond_Urbain(t,x)] = ElCond_Viscosity(T(t),X{:});
        
        %Density
        Comp.Density.Keen(t,x) = Rho_Keen(T(t),X{:});
        Comp.Density.Mills(t,x) = Rho_Mills(T(t),X{:});
        Comp.Density.Stebbins(t,x) = Rho_Stebbins(T(t),X{:});
        Comp.Density.Xin(t,x) = Rho_Xin(T(t),X{:});
        
    end
end

end


