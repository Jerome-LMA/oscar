   % Example of OSCAR script to define a 4 mirrors cavity
clearvars; close all;
addpath(genpath([pwd filesep '..' filesep 'Classes']));

disp('---------------------------------------------------------------------------')
disp('                  OSCAR V3.30                                     ')
disp('  ')


% Define the grid for the simulation: 128 X 128, 3.5 mm X 3.5 mm
G1 = Grid(128,0.15);

% Define the incoming beam before the input mirror surface 
E_input = E_Field(G1,'w',0.02,'R',-1700);

% Define the 4 mirrors of the cavity
I(1) = Interface(G1,'RoC',2500,'CA',0.12,'T',0.001);
I(2) = Interface(G1,'CA',0.12,'T',10E-6);
I(3) = Interface(G1,'CA',0.12,'T',10E-6);
I(4) = Interface(G1,'RoC',2500,'CA',0.12,'T',0.001);

% Distance between the mirrors in m
d = [250 500 250];

% Define the cavity
OMC = CavityN(I,d,E_input);

%OMC.Check_stability
%OMC = Cavity_scan(OMC);
OMC = Cavity_Resonance_Phase(OMC);
OMC = Calculate_Fields_AC(OMC);
Display_Results(OMC)