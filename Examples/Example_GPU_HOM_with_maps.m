clearvars; close all;
addpath(genpath(strcat(pwd, '\..\Classes')));

disp('---------------------------------------------------------------------------')
disp('                  OSCAR V3.30                                      ')
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

IM = Interface(G1,'RoC',1500,'CA',0.4,'T',0.02);
EM = Interface(G1,'RoC',1700,'CA',0.4,'T',2E-6);

% Load the mirror maps
param_PSD = [0.02 -1.4];
Virtual_map_IM = Do_Virtual_Map(G1,param_PSD);
Virtual_map_EM = Do_Virtual_Map(G1,param_PSD);

% Add with 1 nm RMS on the central part
IM = Add_Map(IM,Virtual_map_IM,'reso',G1.Step,'remove_tilt_focus',0.150,'RMS',1E-9,'verbose',false);
EM = Add_Map(EM,Virtual_map_EM,'reso',G1.Step,'remove_tilt_focus',0.150,'RMS',1E-9,'verbose',false);

% To normalise the maps for example to 1 nm RMS use:
% IM = Add_map(IM,'Map1.txt',1.5E-3,1E-9);

% Use the 2 previous Interfaces and the input beam to defing a cavity 3000
% meter long
C1 = Cavity1(IM,EM,3000,E_input);

% To use the digital integration technique
C1.propagation_mat.Use_DI = true;
C1 = Declare_on_GPU(C1);

% Calculate the resonance length
C1 = resonance_phase(C1);

% Calculate and display the reflected field

C1 = calculate_fields_ac(C1); % accelerated convergence scheme
%C1 = Calculate_fields(C1); % old method

display_results(C1);

%% Look at the power content in higher order modes of the circulating field
%figure(2)
%Expand_HOM(C1.field_circ,12,'basis',[0.04275  -1500],'display','vector');

%% Test cavity scan

% C2 = Cavity_Scan(C1);
% C2 = calculate_fields_ac(C2);
% display_results(C2);



