clearvars; close all;
addpath(genpath([pwd filesep '..' filesep 'Classes']));

disp('---------------------------------------------------------------------------')
disp('                  OSCAR V3.30                                     ')
disp('  ')


% Define the grid for the simulation: 256 X 256, 30 cm X 30 cm
G1 = Grid(256,0.40);

% Define one incoming beam (beam radius: 5 cm) 
E_input = E_Field(G1,'w0',0.05);

% Propagate over 100m
E_2 = Propagate_E(E_input,100);

AoI = 45; % Angle of incidence in degree
I1 = Interface(G1,'CA',0.4,'T',0,'RoC',4E3,'AoI',AoI);
[~,E_3] = Transmit_Reflect_Interface(E_2,I1);

% Propagate another time over 100m
E_4 = Propagate_E(E_3,100);

% Then display the beam parameters:
disp('FFT code result:')
Fit_TEM00_xy(E_4);

% If a map has to be inserted, one can use the function:
% Tilt_interface(), suitable to tilt an interface without a strong
% curvature


%% Do the same thing with the ABCD matrix
Lambda = E_input.Wavelength;
q_start = 1/(- 1i*Lambda/(pi*0.05^2));

Mat_propa =  [1 100;0 1]*[1 0;-2/(4E3 * cos (AoI*pi/180) ) 1]*[1 100;0 1];
q_propa = (Mat_propa(1,1)*q_start + Mat_propa(1,2))/(Mat_propa(2,1)*q_start + Mat_propa(2,2));

q_circ_inv = 1/(q_propa);
RofC = 1/real(q_circ_inv);
Beam_rad = sqrt( 1/(-imag(q_circ_inv)*pi/(Lambda)));

disp('---------------------------------------------------')
disp('ABCD matrix result in x:')
fprintf('beam radius: %g      wavefront RofC: %g \n',Beam_rad,RofC)

Mat_propa =  [1 100;0 1]*[1 0;-2/(4E3 / cos (AoI*pi/180) ) 1]*[1 100;0 1];
q_propa = (Mat_propa(1,1)*q_start + Mat_propa(1,2))/(Mat_propa(2,1)*q_start + Mat_propa(2,2));

q_circ_inv = 1/(q_propa);
RofC = 1/real(q_circ_inv);
Beam_rad = sqrt( 1/(-imag(q_circ_inv)*pi/(Lambda)));

disp('ABCD matrix result in y:')
fprintf('beam radius: %g      wavefront RofC: %g \n',Beam_rad,RofC)
