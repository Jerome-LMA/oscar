function Cout = Cavity_resonance_phase(Cin)
%  Cout = Cavity_resonance_phase(Cin) find the resonance of a cavity
% This procedure can find the resonance of the cavity by finding the
% suitable round trip phase shift to bring the circulating field on resonance.

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


num_iter = Cin.Cavity_phase_param;

Field_total = Normalise_E(Field_Circ1,0);
Phase_adjust =1;
Field_Circ = Field_Circ1;

% Calculate the pseudo eigen mode in the cavity
for q = 1:num_iter
    
    Field_total = Field_total + Field_Circ;
    
    % Do a round trip
    Field_Circ = Propagate_E(Field_Circ,Cin.Propagation_mat);
    Field_Circ = Reflect_Mirror(Field_Circ,Cin.I_end);
    Field_Circ = Propagate_E(Field_Circ,Cin.Propagation_mat)*Phase_adjust;
    Field_Circ = Reflect_Mirror(Field_Circ,Cin.I_input);
    
    Phase_adjust = Phase_adjust * exp(-1i*angle(Calculate_Overlap(Field_Circ,Field_total)));
    %     Calculate_power(Field_total)
end

% Then find the round trip to make the eigen mode on resonnance

Field_before = Field_total;
Field_Circ = Field_total;

Field_Circ = Propagate_E(Field_Circ,Cin.Propagation_mat);
Field_Circ = Reflect_Mirror(Field_Circ,Cin.I_end);
%angle(Calculate_Overlap(Field_Circ,Field_before))
Field_Circ = Propagate_E(Field_Circ,Cin.Propagation_mat);
Field_Circ = Reflect_Mirror(Field_Circ,Cin.I_input);

% Field_before.Fit_TEM00;
% Field_Circ.Fit_TEM00;
% angle(Calculate_Overlap(Field_Circ,Field_before))
Cout.Resonance_phase = exp(-1i* angle(Calculate_Overlap(Field_Circ,Field_before)));

disp(['Found the phase for resonance in cavity ' inputname(1)])

% Found an approximation for the resonanting field in the cavity

Cavity_gain = (1-Cin.I_input.r^2)/ (1 - Cin.I_input.r * Cin.I_end.r)^2;
% do the overlap
Coeff_over = Calculate_Overlap(Field_Circ,Field_Circ1);
Pcirc = abs(Coeff_over).^2 * Cavity_gain; % For 1W input power
Cout.Field_reso_guess = Normalise_E(Field_Circ,Pcirc);

% Suppose that the SB are in antiresonance

if Cout.Field_reso_guess.Nb_Pair_SB
    Cavity_gain_SB = (1-Cin.I_input.r^2)/ (1 + Cin.I_input.r * Cin.I_end.r)^2;
    Pcirc_SB = abs(Coeff_over).^2 * Cavity_gain_SB;
    
    tmp_field_SB = Normalise_E(Field_Circ,Pcirc_SB);
    
    for ii=1:Cout.Field_reso_guess.Nb_Pair_SB
        Cout.Field_reso_guess.SB(ii).Field_lower = tmp_field_SB.Field;
        Cout.Field_reso_guess.SB(ii).Field_upper = tmp_field_SB.Field;
    end 
end

end