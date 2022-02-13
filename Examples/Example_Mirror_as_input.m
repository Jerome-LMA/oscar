% Create a cavity with a thick input mirror
% this cavity example follows the parameters of an Advanced Virgo arm
% cavity.

% We will check the etalon effect

clearvars; close all;
addpath(genpath('Classes'));

disp('---------------------------------------------------------------------------')
disp('                  OSCAR V3.30  - The cavities with thick mirror       ')
disp('  ')

G1 = Grid(128,0.5);

E_input = E_Field(G1,'w',0.0491,'R',-1431);


% Arm cavity of Advanced Virgo
% Define the input mirror
% Define the AR side of the ITM, reflection of the AR coating = 100 ppm
IM_AR = Interface(G1,'RoC',-1420,'CA',0.33,'T',1-100E-6,'L',0,'n1',1,'n2',1.45);
% Define the HR side
IM_HR = Interface(G1,'RoC',1420,'CA',0.33,'T',0.014,'L',0,'n1',1,'n2',1.45);
% Finally define a mirror of thickness 20 cm
IM = Mirror(IM_HR,IM_AR,0.2);

% Define the end mirror, only as an interface
EM_HR = Interface(G1,'Roc',1683,'CA',0.33,'T',75E-6,'L',0);

% Create the cavity
C_NA = Cavity1(IM,EM_HR,3000,E_input);
C_NA.resonance_phase()
% Check the cavity parameters
%C_NA.check_stability()

C_NA.resonance_phase();

C_NA.calculate_fields_ac();
C_NA.display_results

 %C_NA.get_info

%% Now, scan the thickness of the input mirror over one wavelength

disp('---------------------------------------------------------------------------')
disp('              OSCAR V3.30     Simulate the etalon effect    ')
disp('  ')



Nb_points = 20;
Scan_length = 1064E-9;

Vec_length = linspace(0,Scan_length,Nb_points);
Vec_ref_power = zeros(Nb_points,1);

f = waitbar(0,'Initialising...');

for ii = 1:Nb_points
    IM = Mirror(IM_HR,IM_AR,0.2+Vec_length(ii));   % Change the length of the substrate
    IM.RT_inside = 3;                              % Do several round trip inside the mirror
    
    C_NA = Cavity1(IM,EM_HR,3000,E_input);
    C_NA = resonance_phase(C_NA,'verbose',false);
    
    C_NA = calculate_fields_ac(C_NA);   
    Vec_ref_power(ii) = calculate_power(C_NA.Field_ref);
    
    waitbar(ii/Nb_points,f,'Scanning the IM thickness');
    
end
close(f)

figure(4)
plot(Vec_length*1E6,Vec_ref_power)
xlabel('IM Thickness variation [µm]')
ylabel('Reflected power [W]')
