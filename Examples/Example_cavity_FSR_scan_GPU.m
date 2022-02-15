clearvars; close all;
addpath(genpath(strcat(pwd, '\..\Classes')));

disp('---------------------------------------------------------------------------')
disp('                  OSCAR V3.21                                      ')
disp('  ')


% Define the grid for the simulation: 128 X 128, 15 cm X 15 cm
G1 = Grid(128,0.15);

% Define the incoming beam on the input mirror surface (beam radius 2 cm,
% wavefront curvature 2500 m, propagating toward the waist)
E_input = E_Field(G1,'w',0.02,'R',-2500);

% Since version 3.17, now possibility to add SB (but it slows the
% calculation)
%E_input = Add_Sidebands(E_input,5E6,0.1);

% Define the 2 mirrors, RofC = 2500m, 10 cm in diameter, transmission 2%,
% no loss

IM = Interface(G1,'RoC',2500,'CA',0.10,'T',0.02);
EM = Interface(G1,'RoC',2500,'CA',0.10,'T',0.02);

% Misaligned the end mirror by 100 microradian
EM = Add_Tilt(EM,1E-6,'y');

% Use the 2 previous Interfaces and the input beam to defing a cavity 1000
% meter long
C1 = Cavity1(IM,EM,1000,E_input);
% Convert EFields to GPU-Arrays
C1 = Declare_on_GPU(C1);

% Calculate the resonance length

C1 = Cavity_Scan(C1,'use_parallel',true,'With_SB',false);

% Display information about the cavity
Display_scan(C1);


