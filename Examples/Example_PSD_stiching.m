clear variables; close all; 
addpath(genpath([pwd filesep '..' filesep 'Classes']));

% Example created with Maxime Lejean during his PhD 2023-2026

disp('---------------------------------------------------------------------------')
disp('                  OSCAR V3.31                                      ')
disp('  ')

% Define the original high resolution map
GL = Grid(4096,0.4);

map = Do_Virtual_Map(GL,'Standard');
I1 = Interface(GL,'RoC',inf);
I1 =  Add_Map(I1,map,'reso',GL.Step);
%E_Plot(I1)

% Resample the map with a lower resolution
G1 = Grid(512,0.4);
I2 = Interface(G1,'RoC',inf);
I2 =  Add_Map(I2,map,'reso',GL.Step);
figure(1);I_Plot(I2)
Weighted_RMS(I2,'diam',0.4)

G1 = Grid(256,GL.Step*256);
I3 = Interface(G1,'RoC',inf);
I3.surface = I1.surface(1:256,1:256);
figure(2);I_Plot(I3)
Weighted_RMS(I3,'diam',GL.Step*256)

[PSD_1D2,Freq_PSD2] = Plot_PSD(I2,'Diam',0.4,'display',false);
[PSD_1D3,Freq_PSD3] = Plot_PSD(I3,'Diam',GL.Step*256,'display',false);

figure(3)
loglog(Freq_PSD2,PSD_1D2, 'linewidth',2)
hold on
loglog(Freq_PSD3,PSD_1D3, 'linewidth',2)
hold off
ylabel('Power Spectral Density [m^3/m^-1]','FontSize',14);
xlabel('Spatial frequency [1/m]','FontSize',14);
legend(' PSD of the full map',' PSD on a zoom part')
grid on; box on;

