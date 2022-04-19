% Example of OSCAR script to define a 4 mirrors cavity and to calculate the
% eigenmodes of the cavity
clearvars; close all;
addpath(genpath([pwd filesep '..' filesep 'Classes']));

disp('---------------------------------------------------------------------------')
disp('                  OSCAR V3.30                               ')
disp('  ')


% Define the grid for the simulation: 256 X 256, 3.5 mm X 3.5 mm

G1 = Grid(256,0.0035);

% Define the 4 mirrors of the cavity
I(1) = Interface(G1,'CA',0.2,'T',0.02);
I(2) = Interface(G1,'CA',0.2,'T',0.02);
I(3) = Interface(G1,'RoC',2.34,'CA',0.0035,'T',0,'AoI',13);    % 'AoI'is for the angle of incidence in degree
I(4) = Interface(G1,'RoC',2.34,'CA',0.0035,'T',10E-6,'AoI',13);

% A dummy input field
E_input = E_Field(G1,'w0',0.02);

% Distance between the mirrors in m
d = [0.07 0.147 0.294 0.147];

% Define the cavity
OMC = CavityN(I,d,E_input);

% Now a trick, with the actual grid size, the round trip matrix to
% calculate the eigenmodes would be too big ( (256 X 256)^2 points) but the high resolution grid is required for the digital
% integration propagation to work.
% So all the calculations for the round trip matrix will be done on the
% high resolution grid, but the results will be downsampled on the smaller
% grid to be sure that the calculation of the eigenmode is manageable
% !! it will be slooooooooow !!

OMC = Calculate_RT_mat(OMC,'Grid',Grid(64,0.0035)); % Notice how a second smaller grid is given, same length as G1
Display_Cavity_Modes(OMC,'N',40);

% Display the Airy peak of the first 30 (lowest clipping loss) eigenmodes:
% Display_cavity_modes(C1,'N',30,'Airy',1)