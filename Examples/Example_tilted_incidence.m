clearvars; close all;
addpath(genpath([pwd filesep '..' filesep 'Classes']));


disp('---------------------------------------------------------------------------')
disp('                  OSCAR V3.31                                     ')
disp('  ')


% Define the grid for the simulation: 256 X 256, 30 cm X 30 cm
G1 = Grid(512,0.40);

% Define one incoming beam (beam radius: 5 cm) 
E_input = E_Field(G1,'w0',0.05);

% Propagate over 100m
E_2 = Propagate_E(E_input,100);

AoI = 10; % Angle of incidence in degree
RoC = 320;

I1 = Interface(G1,'CA',0.4,'T',0,'RoC',RoC,'AoI',AoI);
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

Mat_propa =  [1 100;0 1]*[1 0;-2/(RoC * cos (AoI*pi/180) ) 1]*[1 100;0 1];
q_propa = (Mat_propa(1,1)*q_start + Mat_propa(1,2))/(Mat_propa(2,1)*q_start + Mat_propa(2,2));

q_circ_inv = 1/(q_propa);
RofC = 1/real(q_circ_inv);
Beam_rad = sqrt( 1/(-imag(q_circ_inv)*pi/(Lambda)));

disp('---------------------------------------------------')
disp('ABCD matrix result in x:')
fprintf('beam radius: %g      wavefront RofC: %g \n',Beam_rad,RofC)

Mat_propa =  [1 100;0 1]*[1 0;-2/(RoC / cos (AoI*pi/180) ) 1]*[1 100;0 1];
q_propa = (Mat_propa(1,1)*q_start + Mat_propa(1,2))/(Mat_propa(2,1)*q_start + Mat_propa(2,2));

q_circ_inv = 1/(q_propa);
RofC = 1/real(q_circ_inv);
Beam_rad = sqrt( 1/(-imag(q_circ_inv)*pi/(Lambda)));

disp('ABCD matrix result in y:')
fprintf('beam radius: %g      wavefront RofC: %g \n',Beam_rad,RofC)

%% Use a more complex formula for the wavefront distortion (not implemented in OSCAR by default)
% Hello, P., and J-Y. Vinet. "Simulation of beam propagation in off-axis optical systems." Journal of optics 27, no. 6 (1996): 265.
%https://iopscience.iop.org/article/10.1088/0150-536X/27/6/005

I2 = Interface(G1,'CA',0.4,'T',0,'RoC',inf); % create a dummy interface
theta = AoI*pi/180;
%RoC = -RoC;

f_2D = RoC * cos(theta) - sqrt( (RoC * cos(theta))^2 - 2*RoC * G1.D2_X * sin(theta) - G1.D2_square  );
g_2D = sqrt( (cos(theta))^2 - 2* G1.D2_X * sin(theta) / RoC  - G1.D2_square / RoC^2  );
s_2D = f_2D - ( G1.D2_X * sin(2*theta) - f_2D * cos(2*theta) ) ./ ( 2* sin(2*theta)*(sin(theta+ G1.D2_X/RoC)) .* g_2D - cos(2*theta)*(1 - 2*g_2D.^2) );
% formule 29 from the article

%I2.surface = -0.5*s_2D; % pass from wavefront to surface
% !! The surface also contains the tilt of the mirror !!
% !! it should be removed depending of the simulations !!
% !! to remove it, use Add_Map()

I2 = Add_Map(I2,s_2D,'reso',G1.Step,'scale',0.5,'remove_tilt',G1.Length);

% figure(4)
% imagesc(I1.surface ./ I2.surface)
% axis square

[~,E_5] = Transmit_Reflect_Interface(E_2,I2);

% Propagate another time over 100m
E_6 = Propagate_E(E_5,100);

% Then display the beam parameters:
disp('FFT code result with realistic formula according to the article:')
Fit_TEM00_xy(E_6);
% large deviation could be noted for angle > 10 deg







