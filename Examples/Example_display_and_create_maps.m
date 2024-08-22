clear variables; close all; 
addpath(genpath([pwd filesep '..' filesep 'Classes']));

% PSD calculations and virtual maps generation updated by Maxime Lejean
% as part of his PhD at LMA 2023-2026

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

%% Display the PSD of the map and do a fit

[PSD_1,xFreq] = Plot_PSD(Iout, 'diam', 0.1, 'display', true,'RMS_display',true); 

% ---- We can fit the PSD in 1D ---- %
% ---- Parameter initial conditions for the fit ---- % 
% [Amplitude,Node1,Node2,...,Law1,Law2,...]
 MP_ini  = [200, 1000, 1 ,2 ]; % 2 Segments, with a corner frequency at 1000 Hz and with different slopes before and after

FitFunctionPSD_log = @(MP,Frequencies) log(FitFunctionPSD(Frequencies,MP));
Fitted_Parameters = lsqcurvefit(FitFunctionPSD_log,MP_ini,xFreq,log(PSD_1));%,LB,UB); % New parameters

%% Create a synthetic map based from a parametrised PSD

% Define a new grid
G2 = Grid(532,0.1);

Generated_map = Do_Virtual_Map(G2,Fitted_Parameters);

I3 = Interface(G2,'RoC',Inf);
I3 = Add_Map(I3,Generated_map,'reso',G2.Step);

figure(4);I_Plot(I3,'diam',0.1)

% Comparison of the RMS:
% Weighted_RMS(Iout,'diam',0.1)
% Weighted_RMS(I3,'diam',0.1)


%% Compared the PSDs of the original map vs the synthetic new one and the fitted PSD

% Calculate the fit of the PSD
Fitted_PSD_1 = FitFunctionPSD(xFreq,Fitted_Parameters);

% Calculate the PSD of the new map
[PSD_2,xFreq2] = Plot_PSD(I3, 'diam', 0.1,'display',false); 

figure(5); hold on; 
set(gca, 'XScale', 'log', 'YScale', 'log','fontsize',14); 
loglog(xFreq,PSD_1,'r','LineWidth',2);
loglog(xFreq,Fitted_PSD_1,'-b','LineWidth',2);
loglog(xFreq2,PSD_2,'g','LineWidth',2);
hold off; 
ylabel('Power Spectral Density [m^3/m^-1]','FontSize',14);
xlabel('Spatial frequency [1/m]','FontSize',14);
legend(' PSD of the initial map',' Fitted PSD',' PSD of the virtual map')
grid on; box on;























