%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Main Script to aggregate all the properties of the different compositions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize

% Clear the workspace
clear all;
clc ;

% Define the Temp vector
T = [1800,1850,1900,1950]+273 ; % in Kelvin can be single point or vector


% Compositions in wt%! 
%                     SiO2,  TiO2,   Al2O3, Cr2O3,   FeO,   MgO,   MnO,   CaO,  Na2O,   K2O,     Li2O,CaF2,ZrO2,B2O3,CrO,NiO,Fe2O3,BaO,SrO
LHS_1=              [44.18,  0.79,   26.24,  0.02,  3.04, 11.22,  0.05, 11.62,  2.30,  0.46,       0,0,0,0,0,0,0,0,0]'; % Don't forget to transponse!
LMS_1=              [42.81,  4.62,   14.31,  0.21,  7.87, 18.89,  0.15,  5.94,  4.92,  0.57,       0,0,0,0,0,0,0,0,0]';
Apollo16_64501 =    [45.42,  0.45,   27.85,  0.08,  4.49,  4.39,  0.06, 16.77,   0.4,  0.09,       0,0,0,0,0,0,0,0,0]';
Luna24_24999 =      [44.61,  0.99,   10.77,  0.42, 20.83, 10.97,  0.28, 10.87,  0.23,  0.02,       0,0,0,0,0,0,0,0,0]';

% Initialize the structures
LHS = {};
LMS = {};
Apollo = {};
Luna = {};

%% Do the calculations
% Calculate the composition trajectories
LHS.Trajectory = Traj(LHS_1);
LMS.Trajectory = Traj(LMS_1);
Apollo.Trajectory = Traj(Apollo16_64501);
Luna.Trajectory = Traj(Luna24_24999);
% Add the properties to the data structure
LHS = PropertiesDuringO2Removal(T,LHS);
LMS = PropertiesDuringO2Removal(T,LMS);
Apollo = PropertiesDuringO2Removal(T,Apollo);
Luna = PropertiesDuringO2Removal(T,Luna);

%% Save the files

save('Prop_LHS.mat','LHS');
save('Prop_LMS.mat','LMS');
save('Prop_Apollo.mat','Apollo');
save('Prop_Luna.mat','Luna');

