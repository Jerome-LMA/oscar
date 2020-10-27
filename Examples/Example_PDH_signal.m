clearvars; close all;
addpath(genpath('Classes'));

disp('---------------------------------------------------------------------------')
disp('                  OSCAR V3.30                                 ')
disp('  ')


% Define the grid for the simulation: 128 X 128, 40 cm X 40 cm
G1 = Grid(128,0.4);

% Define the incoming beam on the input mirror surface (beam radius 4.3 cm,
% wavefront curvature -1034 m, propagating toward the waist)
E_input = E_Field(G1,'w',0.043,'R',-1034);

% Add the sidebands to the field E1, frequency of modulation 6.7 MHz,
% modulating index 0.3
E_input = Add_Sidebands(E_input,6.7234E6,0.3);

% Define the 2 mirrors, RofC_IM = 1500m, RofC_IM = 1700m, 30 cm in
% diameter, transmission 5% for the input mirror and 0.5% for the end
% mirror.

IM = Interface(G1,'RoC',1500,'CA',0.35,'T',0.05);
EM = Interface(G1,'RoC',1700,'CA',0.35,'T',0.005);

% Use the 2 previous Interfaces and the input beam to defing a cavity 3000
% meter long
C1 = Cavity1(IM,EM,3000,E_input);
C1.Laser_start_on_input = false ;
%C1.Propagation_mat.Use_DI = true;
C1 = Cavity_Resonance_Phase(C1);
C1 = Calculate_Fields_AC(C1);

Nb_point = 200;                                     % Number of points for the scan

Phase_scan = zeros(Nb_point,1);           % Phase shift used to scan the cavity
Sig.p = zeros(Nb_point,1);                      % Demodulated signal in phase in reflection from the cavity
Sig.q = zeros(Nb_point,1);                      % Demodulated signal in quadrature in reflection from the cavity
Power.car = zeros(Nb_point,1);               % Circulating power of the carrier in the cavity
Power.SBl = zeros(Nb_point,1);               %  Circulating power of the lower sideband in the cavity
Power.SBu = zeros(Nb_point,1);              %  Circulating power of the upper sideband in the cavity

for ii=1:Nb_point
    Pct = round(ii*100/Nb_point);
    fprintf(1,'%3i %%',Pct);
    
    Phase_scan(ii) = ii*(2*pi)/Nb_point;          % Scan the round trip phase shift from 0 to 2 pi
    C1.Resonance_phase = exp(1i*Phase_scan(ii));  % Set the round trip phase shift for the cavity
    
    tic
    C1 = Calculate_Fields_AC(C1);
    C1.Field_reso_guess = C1.Field_circ;
    [Sig.p(ii),Sig.q(ii)] = Demodulate_SB(C1.Field_ref,'phase',pi/2);       % Demodule the carrier with the sidebands in reflection
    time_need(ii) = toc;
    Power.car(ii) = Calculate_Power(C1.Field_circ);           % calculate also the power of the carrier circulating in the cavity
    [Power.SB1(ii),Power.SB2(ii)] = Calculate_Power(C1.Field_circ,'SB');   % and the sidebands
    
    if ii ~= Nb_point
        fprintf(1,'\b\b\b\b\b \b \b \b')
    else
        fprintf(1,'done \n')
    end
end

% Plot all the results
figure(3)
hold all
plot(Phase_scan,Sig.p,'LineWidth',2)
plot(Phase_scan,Sig.q,'LineWidth',2)
grid on; box on;
hold off
legend('Signal in phase','Signal in quadrature')
title('Demodulated PDH signal in reflection')
xlabel('Cavity round trip phase shift')
ylabel('Signal [a.u.]')

%
figure(4)
semilogy(Phase_scan,Power.car,'LineWidth',2)
grid on; box on;
hold all
semilogy(Phase_scan,Power.SB1,'LineWidth',2)
semilogy(Phase_scan,Power.SB2,'LineWidth',2)
hold off
legend('Carrier','Lower sideband','Upper sideband')
title('Power of the fields circulating inside the cavity')
xlabel('Cavity round trip phase shift')
ylabel('Power [W]')
