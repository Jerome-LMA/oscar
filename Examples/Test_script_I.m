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

Dummy = Add_Map(Dummy,Map_loaded,'reso',dx,'remove_tilt_focus',0.200,'shift',[0.0,-0.02]);
figure(2)
I_Plot(Dummy,'diam',0.20)