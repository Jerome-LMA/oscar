clear variables; close all; 
addpath(genpath([pwd filesep '..' filesep 'Classes']));

disp('---------------------------------------------------------------------------')
disp('                  OSCAR V3.31                                      ')
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
%Dummy = Add_Map(Dummy,Map_loaded,'reso',dx,'remove_tilt_focus',0.250,'shift',[-0.02,0.01]);

% Plot the map over a diameter of 25cm
figure(2)
I_Plot(Dummy,'diam',0.25)

% !! very important: this map is equivalent to wavefront distortions as
% seens from the incoming laser beam. It is the convention in OSCAR
% the height is opposite from the surface error, where a bump has a
% positive height.

%% Load directly a map and create the matching grid at the same time 
% fastest way when one has just to plot the map without any simulation, without any interpolation.
% the 2D plot will be similar as seen in the Zygo Metropro interface.
% Sign flip since OSCAR 3.20 compared to previous versions for the function Load_dat_map 
% (change after I noticed how people are using OSCAR)

[Iout,Gout] = Load_dat_Map('Example_ZYGO_data.dat','remove_tilt_focus',0.1);
figure(3);I_Plot(Iout,'diam',0.1)
