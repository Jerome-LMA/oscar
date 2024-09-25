clearvars; close all;
addpath(genpath([pwd filesep '..' filesep 'Classes']));

disp('---------------------------------------------------------------------------')
disp('                  OSCAR V3.31                                     ')
disp('  ')


% Define the grid for the simulation: 256 X 256, 15 cm X 15 cm
G1 = Grid(256,0.15);

% Define the incoming beam outside the input mirror (beam radius 2 cm), at
% the waist

E_input = E_Field(G1,'w0',0.022);
%E_input = Add_Sidebands(E_input,'Mod_freq',3.4E6,'Mod_depth',0.2);

% Imperfect mode matching and the beam will be also be slightly clipped on
% the end mirror (diffraction loss of 1700 ppm)

% Define the 2 mirrors, one flat and the other with a RoC of 2400m, 10 cm in diameter, transmission 2% and 0.1%,
% no loss

IM = Interface(G1,'RoC',inf,'CA',0.1,'T',0.02);
EM = Interface(G1,'RoC',2400,'CA',0.1,'T',0.001);


% Use the 2 previous Interfaces and the input beam to defing a cavity 1000
% meter long
C1 = Cavity1(IM,EM,1000,E_input);

% Calculate the resonance length
C1 = Cavity_Resonance_Phase(C1);

% Display the circulating power, reflected and transmitted powers

tic
C2 = Calculate_Fields_AC(C1);
disp('Accelerated convergence results:')
AC_time = toc;
C2.Display_Results('display',false);


tic
C3 = Calculate_Fields(C1,'accuracy',0.00001);
disp('Normal convergence results:')
NC_time = toc;
C3.Display_Results('display',false);

fprintf('\n Computational speed gain: %3.2g \n', NC_time/AC_time )

