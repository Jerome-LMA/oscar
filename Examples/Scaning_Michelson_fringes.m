% Script as a pedagogical example for a scanning Michelson interferometer
% with a flat mirror (reference) in one arm and a curved mirror into the other arm
% (whom we could measure the radius of curvature)

clearvars; close all;
addpath(genpath(strcat(pwd, '\..\Classes')));

disp('---------------------------------------------------------------------------')
disp('                  OSCAR V3.20                                     ')
disp('  ')


G1 = Grid(256,0.10);

% Define the input beam
E_input = E_Field(G1,'Wavelength',633E-9,'w0',0.04);

E1 = Propagate_E(E_input,0.25);

Ref_BS = 0.5; % reflection of the beamsplitter in power

% Split between 2 arms
E_north = sqrt(Ref_BS) * Propagate_E(E1,0.25);
E_east = 1i*sqrt(1-Ref_BS) * Propagate_E(E1,0.25);

% Propagate the light between the 2 arms
E_north = Propagate_E(E_north,0.25);
E_east = Propagate_E(E_east,0.25);

% Reflect, let's say the north arm has a flat reference mirror
E_north  = Reflect_Mirror(E_north,1E99);
E_east  = Reflect_Mirror(E_east,600);

E_north = 1i*sqrt(1-Ref_BS) * Propagate_E(E_north,0.25);
E_east = sqrt(Ref_BS) * Propagate_E(E_east,0.25);

Vec_phase = linspace(0,4*pi,100);

for ii = 1:length(Vec_phase)
    % Propagate back and then interfere   
    Eout = E_north*exp(1i*Vec_phase(ii)) + E_east;
    
    E_Plot(Eout,'display','intensity');caxis([0 300]);colormap(hot);pause(0.05);
end


