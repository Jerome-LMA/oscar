clear variables; close all; clear classes
addpath(genpath(strcat(pwd, '\..\Classes')));

disp('---------------------------------------------------------------------------')
disp('                  OSCAR V3.30                                      ')
disp('  ')

%Define the grid for the simulation: 800 X 800, 40 cm X 40 cm
G1 = Grid(800,0.4);

%--------------------------------------------------------------------------------------------
% Tips to check the look of the map loaded

% Create a dummy flat mirror
Dummy = Interface(G1,'RoC',Inf);

% Load the mirror map as a ZYGO.dat file
[Map_loaded, dx] = ReadZygoBinary('Example_ZYGO_data.dat');
figure(1)
imagesc(Map_loaded); axis square % the raw map

% Add the loaded map to the flat mirror, could be used in the simulations
% of a cavity.
Dummy = Add_Map(Dummy,Map_loaded,'reso',dx,'remove_tilt_focus',0.250);

% Plot the map over a diameter of 25cm
figure(2)
I_Plot(Dummy,'diam',0.25)

% !! very important: this map is equivalent to wavefront distortions as
% seens from the incoming laser beam. It is the convention in OSCAR
% the height is opposite from the surface error, where a bump has a
% positive height.
% To plot the map with the physical height, look in the next section

%% Load directly a map and create the matching grid at the same time 
% fastest way when one has just to plot the map without any simulation, without any interpolation.
% the 2D plot will be similar as seen in the Zygo Metropro interface.
% Sign flip in OSCAR 3.20 compared to previous version for the function Load_dat_map 
% (change after I noticed how people are using OSCAR)

[Iout,Gout] = Load_dat_Map('Example_ZYGO_data.dat','remove_tilt_focus',0.1);
figure(3);I_Plot(Iout,'diam',0.1)

%% Create a synthetic map based on a parametrised PSD
% All the credit and merit go to F. Bondu (Virgo) for the function to
% generate the virtual maps from a 1D PSD. I only wrapped the functions to be
% integrated with OSCAR

% Define a grid
G2 = Grid(532,0.4);

% Define the power law for the PSD
% Law was derived according to the various wavefront maps measured at LMA over the years.
% The maps are the ones of large optics (> 300 mm diameter) with nm
% flatness over the central part, so pretty good polishing
% Do not overanalysed the PSD or draw conclusions from them as they
% dependent on a lot of conditions and could change over time. 
% It is there for only illustrative purpose, absolute no warranty that
% you wll have it your optics!
% Do not forget to add which the flatness you are are looking for

param_PSD_ZYGO_IBF = [0.06 2.4 0.004 -4.3 -0.08 -3.2 16 420]; % approximated parameters for ZYGO polishing with IBF (Richmond site), derived from all the AdV IM and EM wavefront measurements
param_PSD_ZYGO_F = [0.08 -2.5]; % approximated parameters for ZYGO polishing (Middlefield site), flat surface, derived from the AdV BS, CP wavefront measurements
param_PSD_ZYGO_C = [3E-4 10 -1 -4 150]; % approximated parameters for ZYGO polishing (Richmond site), curved surface, derived from AdV PR, SR wavefront measurements
param_PSD_Coastline = [0.02 -1.4]; % approximated parameters for Coastline Optics 
param_PSD_General_Optics = [0.03 -2]; % approximated arameters for General_Optics polishing (Gooch & Housego now)

fake_map = Do_Virtual_Map(G2,param_PSD_ZYGO_IBF);

% Added for a flat interface
I1 = Interface(G2,'RoC',Inf);
I1 = Add_Map(I1,fake_map,'reso',G2.Step,'remove_tilt_focus',0.250,'RMS',0.5E-9,'verbose',false);

figure(4);I_Plot(I1,'diam',0.25)
Weighted_RMS(I1,'diam',0.25);




















