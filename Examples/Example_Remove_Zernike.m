clearvars; close all;
addpath(genpath(strcat(pwd, '\..\Classes')));

disp('---------------------------------------------------------------------------')
disp('                  OSCAR V3.20                                      ')
disp('  ')


%Define the grid for the simulation: 256 X 256, 40 cm X 40 cm
G1 = Grid(800,0.4);

%--------------------------------------------------------------------------------------------
% Tips to check the look of the map loaded

% Load the mirror map as a ZYGO.dat file
[Map_loaded, dx] = ReadZygoBinary('example_ZYGO_data.dat');
% Calculate the diameter of the map
map_side = size(Map_loaded)*dx;

% Create a dummy flat mirror
I1 = Interface(G1,'RoC',Inf,'CA',min(map_side));

% Add the loaded map to the flat mirror
I1 = Add_Map(I1,Map_loaded,'reso',dx,'scale',-1);

% Plot the map over a diameter of 25cm
figure(1)
I_Plot(I1,'diam',0.25)

% Alternative way, shorter:
[I2,G2] = Load_dat_Map('example_ZYGO_data.dat');
figure(2); I_Plot(I2,'diam',0.25)

%-----------------------------------------------------------------------------------------
% Remove now the Zernike polynomials over a diameter of 100mm
[I3,test] = Expand_Zernike(I1,'Z_order',12,'diam',0.100);
% Now the reconstructed surface is in the variable I2 (only define over 0.1 m)

% Plot the high frequency: surface -  reconstructed surface up to the order
% 12
figure(3)
I_Plot(I1-I3,'diam',0.099); axis square
caxis([-1E-9,1E-9])

