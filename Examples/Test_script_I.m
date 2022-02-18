clearvars; close all;
addpath(genpath(strcat(pwd, '\..\Classes')));

disp('---------------------------------------------------------------------------')
disp('                  OSCAR V3.21                                  ')
disp('  ')

%% Test Add_Map()
disp('----------   Test Add_Map  -----------')

G1 = Grid(256,0.4);
Dummy = Interface(G1,'RoC',Inf);
[Map_loaded, dx] = ReadZygoBinary('Example_ZYGO_data.dat');
figure(1)
imagesc(Map_loaded); axis square % the raw map

I2 = Add_Map(Dummy,Map_loaded,'reso',dx,'remove_tilt_focus',0.200,'shift',[0.0,-0.02]);
figure(2)
I_Plot(I2,'diam',0.20)

%% Test Add_Astigmatism()

Dummy = Interface(G1,'RoC',Inf);
I3 = Add_Astigmatism(Dummy,2E-9,0.150);
I_Plot(I3,'diam',0.20)

%% Test Cut_frequency_Interface()

I4 = Cut_Frequency_Interface(Dummy,'HP',50);
figure(2)
I_Plot(I4,'diam',0.20)

%% Test Expand_Zernike()

Expand_Zernike(I2,'Z_order',25,'diam',0.2)

I5 = Expand_Zernike(I2,'Z_order',5,'diam',0.2);
I_Plot(I5,'diam',0.20)

%% Test Plot_PSD()

Plot_PSD(I2,'diam',0.2);
[PSD_1D,freq] = Plot_PSD(I2,'diam',0.2,'display',false);

%% Test plus minus
I7 = I5 + I4;
I8 = I3 - I2;

%% 




