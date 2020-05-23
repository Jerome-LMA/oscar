clearvars; close all;
addpath(genpath(strcat(pwd, '\..\Classes')));

disp('---------------------------------------------------------------------------')
disp('                  OSCAR V3.20                                      ')
disp('  ')


% Define the grid for the simulation: 256 X 256, 10 cm X 10 cm
G1 = Grid(256,0.10);

% Define one incoming beam (beam radius: 5 cm) 
E_input = E_Field(G1,'w',0.0125);

% Propagate over 100m
E_2 = Propagate_E(E_input,0.05);

% Then define a telescope as:
% 1 first lens of focal length 8.2 m
% Propagate 4.45 m
% 1 second lens of focal length -3.6 m
% Then propagate 0.6 m

[E3,G3] = Focus_beam_with_telescope(E_2,[8.2 4.45 -3.6 0.6]);

% Then display the beam parameters:
disp('FFT code result:')
Fit_TEM00(E3)
% One can also look at the new lenght of the grid
G3.Length