% Example of OSCAR script to define a 4 mirrors cavity
clearvars; close all;
addpath(genpath([pwd filesep '..' filesep 'Classes']));

disp('---------------------------------------------------------------------------')
disp('                  OSCAR V3.30                                    ')
disp('  ')


% Define the grid for the simulation: 128 X 128, 3.5 mm X 3.5 mm
G1 = Grid(128,0.0035);

% Define the incoming beam before the input mirror surface (beam waist 0.450 mm,
% we start 35 mm before the waist)
E_input = E_Field(G1,'w0',450E-6,'Z',-0.035);

% Define the 4 mirrors of the cavity
I(1) = Interface(G1,'CA',0.2,'T',0.02);
I(2) = Interface(G1,'CA',0.2,'T',0.02);
I(3) = Interface(G1,'RoC',2.34,'CA',0.02,'T',10E-6,'AoI',13);    % 'AoI'is for the angle of incidence in degree
I(4) = Interface(G1,'RoC',2.34,'CA',0.02,'T',10E-6,'AoI',13);

% Distance between the mirrors in m
d = [0.07 0.147 0.294 0.147];

% Define the cavity
OMC = CavityN(I,d,E_input);

%OMC = Cavity_Scan(OMC,'use_parallel',false,'Define_L_length',true);
%Display_Scan(OMC,'scan','RT phase');

% Calculate the resonance
OMC = Cavity_Resonance_Phase(OMC);

% Calculate the steady state fields
OMC = Calculate_Fields_AC(OMC); % Fast convergence scheme (not not implemented with the new SB scheme))
%OMC = Calculate_fields(OMC); % The more traditional method (slower but can be used with sidebands)
Display_Results(OMC)

