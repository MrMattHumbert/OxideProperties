%{
This function maps the trajectory of the 
%}

function traj = Traj(comp_mass)

% Evolution of g O2 removed

%Convert Wt% to Molfraction
%             SiO2, TiO2, Al2O3, Cr2O3,   FeO,   MgO, MnO,   CaO, Na2O,  K2O, Li2O, CaF2, ZrO2, B2O3, CrO, NiO, Fe2O3,   BaO,   SrO
mol_weight = [  60, 79.9, 101.9,   152, 71.85,  40.3,  71, 56.08,   62, 94.2, 29.8,   78,123.2, 69.8,  68,74.7, 159.7, 153.3, 103.6]';
comp_mass = comp_mass./(sum(comp_mass)/100); %Normalize in case something is messed up. 
if (sum(comp_mass)-100) >= 1
    error('The sum of weight fractions does not equal 100');
end
comp_mol = (comp_mass./mol_weight);

maxO2_mol = comp_mol(1)+ comp_mol(4)*2/3 + comp_mol(5)/2 +comp_mol(7)/2 ; %Add the O2 from the SiO2,Cr2O3,FeO,MnO
maxO2_g = maxO2_mol*32; %MW O2 is 32g/mol

O2_Removed = linspace(0,maxO2_g,10000);
del_x = O2_Removed(2)-O2_Removed(1);
comp_mass_evo(:,1) = comp_mass;
for i=2:length(O2_Removed)
    comp_mass_evo(:,i)=comp_mass_evo(:,i-1); %set the current composition to the previous one prior to manipulation
    if comp_mass_evo(5,i-1)>0 % if there is iron remove some
        comp_mass_evo(5,i)=comp_mass_evo(5,i-1)-del_x*mol_weight(5)/32;
    elseif comp_mass_evo(4,i-1)>0 && comp_mass_evo(5,i)<=0 %Remove Chromia
        comp_mass_evo(5,i)=0; % Set iron oxide to zero 
        comp_mass_evo(4,i)=comp_mass_evo(4,i-1)-del_x*mol_weight(4)/32;
    elseif comp_mass_evo(7,i-1)>0 && comp_mass_evo(5,i)<=0 && comp_mass_evo(4,i)<=0 %Remove Manganese
        comp_mass_evo(5,i)=0; % Set iron oxide to zero 
        comp_mass_evo(4,i)=0; % Set Chromia to zero
        comp_mass_evo(7,i)=comp_mass_evo(7,i-1)-del_x*mol_weight(7)/32;
    elseif comp_mass_evo(1,i-1)>0 && comp_mass_evo(7,i)<=0 && comp_mass_evo(5,i)<=0 && comp_mass_evo(4,i)<=0 % set all to zero and remove some silicon
        comp_mass_evo(5,i)=0; % Set iron oxide to zero 
        comp_mass_evo(4,i)=0; % Set Chromia to zero
        comp_mass_evo(7,i)=0; % Set Manganese to zero
        comp_mass_evo(1,i)=comp_mass_evo(1,i-1)-del_x*mol_weight(1)/32; 
    elseif comp_mass_evo(1,i-1)<0
        comp_mass_evo(1,i)=0;
    end
end

% Populate the Structure with the composition Trajectory
traj.MaxO2_g = maxO2_g;
traj.MaxO2_mol = maxO2_mol;
traj.MassEvolution = comp_mass_evo;
traj.MolEvolution = comp_mass_evo./mol_weight;
traj.MolFracEvolution = traj.MolEvolution./sum(traj.MolEvolution,1);

end

