clearvars; close all;
addpath(genpath([pwd filesep '..' filesep 'Classes']));

disp('---------------------------------------------------------------------------')
disp('                  OSCAR V3.31                                     ')
disp('  ')


% Define the grid for the simulation: 256 X 256, 10 cm X 10 cm
G1 = Grid(256,0.10);

% Define one incoming beam (beam radius: 5 cm) 
E_input = E_Field(G1,'w',0.0125);

% Propagate over 10m
E2 = Propagate_E(E_input,10);

% Then define a telescope as:
% 1 first lens of focal length 8.2 m
% Propagate 4.45 m
% 1 second lens of focal length -3.6 m
% Then propagate 0.6 m

[E3,G3] = Focus_Beam_With_Telescope(E2,[8.2 4.45 -3.6 0.6]);

% Then display the beam parameters:
disp('FFT code result:')
Fit_TEM00(E3)
% One can also look at the new lenght of the grid
fprintf('New length of the grid: %g \n',G3.Length)

%% Check with the ABCD matrix

Lambda = 1064E-9;
q_start = 1/(- 1i*Lambda/(pi*0.0125^2));

% Propagate over 10m
Mat_propa = [1 10;0 1];

% Then define a telescope as:
% 1 first lens of focal length 8.2 m
% Propagate 4.45 m
% 1 second lens of focal length -3.6 m
% Then propagate 0.6 m

Mat_propa = [1 0.6;0 1]*[1 0;1/3.6 1]*[1 4.45;0 1]*[1 0;-1/8.2 1]*Mat_propa;
q_propa = (Mat_propa(1,1)*q_start + Mat_propa(1,2))/(Mat_propa(2,1)*q_start + Mat_propa(2,2));

 q_circ_inv = 1/(q_propa);
 RofC = 1/real(q_circ_inv);
 Beam_rad = sqrt( 1/(-imag(q_circ_inv)*pi/(Lambda)));
 
 disp('ABCD matrix result:') 
 fprintf('beam radius: %g      wavefront RofC: %g \n',Beam_rad,RofC)









