clearvars; close all;
addpath(genpath(strcat(pwd, '\..\Classes')));

disp('---------------------------------------------------------------------------')
disp('                  OSCAR V3.20                                   ')
disp('  ')


% Define the grid for the simulation: 128 X 128, 30 cm X 30 cm
G1 = Grid(256,0.1);

% Define the incoming beam (beam radius 2 cm, at the waist)
E1 = E_Field(G1,'w0',0.02);

% Define the biconvex convergent thick lens made of fused silica (n=1.45),
% the length of the substrate is 4 cm

% First surface RofC = 1000 m, 10 cm of aperture, transmission = 1
L1 = Interface(G1,'RoC',-1000,'CA',0.10,'T',1,'n1',1,'n2',1.45);
% Second surface RofC = 2000 m, 8 cm of aperture, transmission = 1
L2 = Interface(G1,'RoC',2000,'CA',0.10,'T',1,'n1',1.45,'n2',1);

E2 = Transmit_Reflect_Interface(E1,L1);
E2 = Propagate_E(E2,0.04);
E2 = Transmit_Reflect_Interface(E2,L2);

E3 = Propagate_E(E2,200);

% Display the beam parameter
Fit_TEM00(E3);

%% Now using the mirror class

% Define the first surface arbitrary called HR, then the second AR
Thick_L = Mirror(L1,L2,0.04);

% Does not really matter the direction of the lens, we can start with the
% HR or AR side.
E2 = Transmit_Reflect_Mirror(E1,Thick_L,'HR');

E3 = Propagate_E(E2,200);

% Display the beam parameter
disp('With the mirror class:') 
Fit_TEM00(E3);


%% Compare the results to ABCD matrix
 
% Define the input laser beam

Lambda = 1064E-9;
q_start = 1/(- 1i*Lambda/(pi*0.02^2));

Mat_propa = [1 200;0 1]*[1 0;((1-1.45)/1)*1/2000 1.45]*[1 0.04;0 1]*[1 0;-((1.45-1)/1.45)*1/1000 1/1.45];

 q_propa = (Mat_propa(1,1)*q_start + Mat_propa(1,2))/(Mat_propa(2,1)*q_start + Mat_propa(2,2));
 
  q_circ_inv = 1/(q_propa);
 RofC = 1/real(q_circ_inv);
 Beam_rad = sqrt( 1/(-imag(q_circ_inv)*pi/(Lambda)));
 
 disp('ABCD matrix result:') 
 fprintf('beam radius: %g      wavefront RofC: %g \n',Beam_rad,RofC)