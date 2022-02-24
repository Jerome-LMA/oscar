clearvars; close all;
addpath(genpath([pwd filesep '..' filesep 'Classes']));

disp('---------------------------------------------------------------------------')
disp('                  OSCAR V3.30                                              ')
disp('  ')


% Define the grid for the simulation: 80 X 80, 35 cm X 35 cm
G1 = Grid(80,0.35);

% Define the incoming beam, it is a dummy beam since it is not used to
% calculate the cavity eigen modes. However it is necessary to define the
% cavity
E_input = E_Field(G1,'w',0.02);

% Advanced Virgo
IM = Interface(G1,'RoC',1420,'CA',0.33,'T',0.01);
EM = Interface(G1,'RoC',1672,'CA',0.33,'T',0.01);
C1 = Cavity1(IM,EM,3000,E_input);

% Advanced LIGO
% IM = Interface(G1,'RoC',1934,'CA',0.33,'T',0.01);
% EM = Interface(G1,'RoC',2245,'CA',0.33,'T',0.01);
% C1 = Cavity1(IM,EM,4000,E_input);

% Calculate the kernel of the cavity
C1 = Calculate_RT_mat(C1);

% Display the cavity eigen modes
Display_Cavity_Modes(C1,'N',30);
  
%Display the Airy peak of the first 30 (lowest clipping loss) eigenmodes with the list of the modes:
%Display_cavity_modes(C1,'N',10,'Airy',true,'List',true)

% Display a particular eigen mode:
% Tab_EM = Display_cavity_modes(C1,'N',20);
% figure(2); E_plot(Tab_EM(5));
