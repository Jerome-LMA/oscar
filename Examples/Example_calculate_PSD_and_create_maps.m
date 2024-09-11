clear variables; close all; 
addpath(genpath([pwd filesep '..' filesep 'Classes']));

% Example created with Maxime Lejean during his PhD 2023-2026

disp('---------------------------------------------------------------------------')
disp('                  OSCAR V3.31                                      ')
disp('  ')

%% Display the PSD of the map and do a fit
% first we have to load a map

[I1,G1] = Load_dat_Map('Example_ZYGO_data.dat','remove_tilt_focus',0.1);
figure(1);I_Plot(I1,'diam',0.1)

[PSD_1,xFreq] = Plot_PSD(I1, 'diam', 0.1, 'display', true,'RMS_display',false); 

% ---- We can fit the PSD in 1D ---- %
% ---- Parameter initial conditions for the fit ---- % 
% [Amplitude,Node1,Node2,...,Law1,Law2,...]
 MP_ini  = [200, 1000, 2 ,4 ]; % 2 Segments, with a corner frequency at 1000 Hz and with different slopes before and after

FitFunctionPSD_log = @(MP,Frequencies) log(FitFunctionPSD(Frequencies,MP));
Fitted_Parameters = lsqcurvefit(FitFunctionPSD_log,MP_ini,xFreq,log(PSD_1));%,LB,UB); % New parameters

%% Create a synthetic map based from a parametrised PSD

% Define a new grid, a larger one to test the extrapolation
G2 = Grid(512,0.2);

Generated_map = Do_Virtual_Map(G2,Fitted_Parameters);

I3 = Interface(G2,'RoC',Inf);
I3 = Add_Map(I3,Generated_map,'reso',G2.Step);

figure(4);I_Plot(I3,'diam',0.2)

%% Compared the PSDs of the original map vs the synthetic new one and the fitted PSD

% Calculate the fit of the PSD
Fitted_PSD_1 = FitFunctionPSD(xFreq,Fitted_Parameters);

% Calculate the PSD of the new map
[PSD_2,xFreq2] = Plot_PSD(I3, 'diam', 0.2,'display',false); 

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























