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
E_input = E_Field(G1,'w0',5E-2,'Include_Birefringence',true);


