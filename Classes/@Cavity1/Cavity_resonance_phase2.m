function Cout = Cavity_resonance_phase2(Cin)
%  Cout = Cavity_resonance_phase(Cin) find the resonance of a cavity
% This procedure can find the resonance of the cavity by finding the
% suitable round trip phase shift to bring the circulating field on resonance.
% More accurate than Cavity_resonance_phase(Cin)

Cout = Cin;

% Normalise the input to 1W
Cin.Laser_in = Normalise_E(Cin.Laser_in);

if ~Cin.Laser_start_on_input
    if isa(Cin.I_input, 'Interface')
        Cin.Laser_in =  Change_E_n(Cin.Laser_in,Cin.I_input.n2);
    end
    Field_Circ1 = Transmit_Reflect_Optic(Cin.Laser_in,Cin.I_input);
else
    Field_Circ1 = Cin.Laser_in;
end

FT_iter = 40;      % Number of iteration for the fine tuning
RT_iter = 40;      % Number of round trip to do for each iteration

% Do a first iteration to have an idea of the resonance length for the
% input beam (supposed to be not far from the expected circulating)

Field_circ_mean = Normalise_E(Field_Circ1,0);
Field_circ_after = Field_Circ1;

Nb_RT = 100;

% for pp = 1:Nb_RT;
%     %figure(1);E_plot(Field_circ_after);pause(0.1)
%
%     Field_circ_before = Field_circ_after;
%     Field_Circ = Propagate_E(Field_circ_before,Cin.Propagation_mat);
%     Field_Circ = Reflect_mirror(Field_Circ,Cin.I_end);
%     Field_Circ = Propagate_E(Field_Circ,Cin.Propagation_mat);
%     Field_circ_after = Reflect_mirror(Field_Circ,Cin.I_input);
%
%     Field_circ_mean = Field_circ_mean + Field_circ_after;
%
%     angle_reso_vec(pp) = angle(Calculate_Overlap(Field_circ_before,Field_circ_after));
% end
Field_total = Normalise_E(Field_Circ1,0);
Field_Circ2 = Field_Circ1;
Phase_adjust = 1;

for q = 1:Nb_RT
    
    %E_plot(Field_total); pause()
    Field_total = Field_total + Field_Circ2;
    
    % Do a round trip
    Field_Circ = Propagate_E(Field_Circ2,Cin.Propagation_mat);
    Field_Circ = Reflect_mirror(Field_Circ,Cin.I_end);
    Field_Circ = Propagate_E(Field_Circ,Cin.Propagation_mat)*Phase_adjust;
    Field_Circ = Reflect_mirror(Field_Circ,Cin.I_input);
    
    %  tmp_vec(q) = angle(Calculate_Overlap(Field_Circ2,Field_Circ));
    Field_Circ2 = Field_Circ;
    Phase_adjust = Phase_adjust * exp(-1i*angle(Calculate_Overlap(Field_Circ,Field_total)));
end

Field_Circ = Propagate_E(Field_total,Cin.Propagation_mat);
Field_Circ = Reflect_mirror(Field_Circ,Cin.I_end);
Field_Circ = Propagate_E(Field_Circ,Cin.Propagation_mat);
Field_Circ = Reflect_mirror(Field_Circ,Cin.I_input);

angle_reso = angle(Calculate_Overlap(Field_total,Field_Circ));


% figure(1)
% plot(tmp_vec); pause;
%Field_circ_mean = 0.1 * Field_circ_mean;
%angle_reso = mean(tmp_vec)

% figure(2)
% E_plot(Field_circ_mean); pause;

%Field_circ_before = Field_Circ1;

Field_Circ = Normalise_E(Field_total,Calculate_power(Field_Circ1));



for pp = 1:FT_iter
    
    %Field_Circ = Field_Circ1;
    %Field_Circ = Field_circ_after;
    %Field_Circ = (1/Nb_RT)*Field_circ_mean;
    
    for nn = 1: RT_iter
        %figure(1);E_plot(Field_Circ); title(num2str(nn));pause(0.1)
        %Calculate_power(Field_Circ)
        Field_Circ = Propagate_E(Field_Circ,Cin.Propagation_mat);
        Field_Circ = Reflect_mirror(Field_Circ,Cin.I_end);
        Field_Circ = Propagate_E(Field_Circ,Cin.Propagation_mat)*exp(1i*angle_reso);
        Field_Circ = Reflect_mirror(Field_Circ,Cin.I_input);
        Field_Circ = Field_Circ + Field_Circ1;
    end
    
    Field_circ_before = Field_Circ;
    Field_Circ2 = Propagate_E(Field_circ_before,Cin.Propagation_mat);
    Field_Circ2 = Reflect_mirror(Field_Circ2,Cin.I_end);
    Field_Circ2 = Propagate_E(Field_Circ2,Cin.Propagation_mat);
    Field_circ_after = Reflect_mirror(Field_Circ2,Cin.I_input);
    
    angle_reso = angle(Calculate_Overlap(Field_circ_before,Field_circ_after));
    tmp_vec(pp) = angle_reso;
    
    %figure(1);E_plot(Field_circ_after); title(num2str(pp));pause(1)
    
end

% figure(1)
% plot(tmp_vec); pause;

Cout.Resonance_phase = exp(1i* angle_reso);

disp(['Found the phase for resonance in cavity ' inputname(1)])

Loss_RTL = (Calculate_power(Field_circ_before) - Calculate_power(Field_circ_after))/Calculate_power(Field_circ_before);


% Found an approximation for the resonanting field in the cavity

Cavity_gain = (-Cin.I_input.t^2)/ (1 - sqrt(1-Loss_RTL))^2; % check here
% do the overlap
Coeff_over = Calculate_Overlap(Field_circ_after,Field_Circ1);
Pcirc = abs(Coeff_over).^2 * Cavity_gain; % For 1W input power
%abs(Coeff_over).^2;

% Add the proper phase shift
Phase_shift = angle(Calculate_Overlap(Field_Circ1,Field_Circ));
Field_Circ = Field_Circ * exp(1i*Phase_shift);

Cout.Field_reso_guess = Normalise_E(Field_circ_after,Pcirc);



% Suppose that the SB are in antiresonance

if ~isempty(Cin.Laser_in.Field_SBl)
    
    Cavity_gain_SB = (1-Cin.I_input.r^2)/ (1 + Cin.I_input.r * Cin.I_end.r)^2;
    Pcirc_SB = abs(Coeff_over).^2 * Cavity_gain_SB;
    
    tmp_field_SB = Normalise_E(Field_Circ,Pcirc_SB);
    
    Cout.Field_reso_guess.Field_SBl = tmp_field_SB.Field;
    Cout.Field_reso_guess.Field_SBu = tmp_field_SB.Field;
    
    Cout.Field_reso_guess.Frequency_Offset = Cin.Laser_in.Frequency_Offset;
    
end





