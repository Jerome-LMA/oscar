clearvars; close all;
addpath(genpath(strcat(pwd, '\..\Classes')));

disp('---------------------------------------------------------------------------')
disp('                  OSCAR V3.20                                      ')
disp('  ')


% Define the grid for the simulation: 256 X 256, 40 cm X 40 cm
G1 = Grid(256,0.4);

% Define the incoming beam outside the cavity (beam radius 4.3 cm,
% wavefront curvature 1034 m, mode Laguerre Gauss 3,3)
E_input = E_Field(G1,'w',0.043,'R',-1034,'mode','LG 3 3');

% Define the 2 mirrors, RofC_IM = 1500m, RofC_IM = 1700m, 30 cm in
% diameter, transmission 2% for the input mirror, almost perfectly
% reflective for the end 
% no loss

IM = Interface(G1,'RoC',1500,'CA',0.35,'T',0.02);
EM = Interface(G1,'RoC',1700,'CA',0.35,'T',2E-6);

% Load the mirror maps
IM = Add_map(IM,'Map1.txt','reso',2E-3);
EM = Add_map(EM,'Map2.txt','reso',2E-3);

% To normalise the maps for example to 1 nm RMS use:
% IM = Add_map(IM,'Map1.txt',1.5E-3,1E-9);

% Use the 2 previous Interfaces and the input beam to defing a cavity 3000
% meter long
C1 = Cavity1(IM,EM,3000,E_input);
C1.Laser_start_on_input = false ;

% To use the digital integration technique
C1.Propagation_mat.Use_DI = true;

% Calculate the resonance length
C1 = Cavity_resonance_phase(C1);

% Calculate and display the reflected field

C1 = Calculate_fields_AC(C1); % accelerated convergence scheme
%C1 = Calculate_fields(C1); % old method

Display_results(C1);

%% Look at the power content in higher order modes of the circulating field
figure(2)
Expand_HOM(C1.Field_circ,12,'basis',[0.04275  -1500],'display','vector');






