clearvars; close all;
addpath(genpath('Classes'));

disp('---------------------------------------------------------------------------')
disp('                  OSCAR V3.30                                   ')
disp('  ')


% Define the grid for the simulation: 256 X 256, 40 cm X 40 cm
G1 = Grid(256,0.4);

% Example arm cavity Advanced Virgo

% Define the input beam, instead of giving some beam parameters, just
% assume perfect mode matching wih the cavity
E_input = E_Field(G1,'Optimal_MM',true);
E_input = Add_Sidebands(E_input,'Mod_freq',10E6,'Mod_index',0.1);

% Define the input and end mirrors of the cavity
IM = Interface(G1,'RoC',1420,'CA',0.33,'T',0.014);
EM = Interface(G1,'RoC',1683,'CA',0.33,'T',5E-6);

% Use the 2 previous Interfaces and the input beam to defing a cavity 3000
% meter long
Arm_cavity = Cavity1(IM,EM,3000,E_input);

% Calculate the resonance length to maximise the circulating power
Arm_cavity = Cavity_Resonance_Phase(Arm_cavity);

% Display information about the cavity
Arm_cavity = Calculate_Fields_AC(Arm_cavity);

Arm_cavity.Display_Results


