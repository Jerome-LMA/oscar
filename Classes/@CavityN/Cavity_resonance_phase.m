function Cout = Cavity_Resonance_Phase(Cin)
%  Cout = Cavity_resonance_phase(Cin) find the resonance of a cavity
% This procedure can find the resonance of the cavity by finding the
% suitable round trip phase shift to bring the circulating field on resonance.

Cout = Cin;

% Normalise the input to 1W
Cin.Laser_in = Normalise_E(Cin.Laser_in);

if ~Cin.Laser_start_on_input
    Field_in =  Change_E_n(Cin.Laser_in,Cin.I_array(1).n2);
    Field_in = Transmit_Reflect_Interface(Field_in,Cin.I_array(1));
    Field_Circ = Field_in;
else
    Field_Circ = Cin.Laser_in;
end

Field_in2 = Field_Circ;

num_iter = Cin.Cavity_phase_param;

Field_total = Normalise_E(Field_Circ,'Power',0);
Phase_adjust = 1;


% Calculate the pseudo eigen mode in the cavity

for q = 1:num_iter
    
    Field_total = Field_total + Field_Circ;
    
    if Cin.type == 'ring'
        for pp=1:Cin.Nb_mirror
            
            if pp ~= Cin.Nb_mirror % check we are not at the last iteration
                Field_Circ = Propagate_E(Field_Circ,Cin.Propagation_mat_array(pp));
                Field_Circ = Reflect_Mirror(Field_Circ,Cin.I_array(pp+1));
            else % we are at the last iteration
                Field_Circ = Propagate_E(Field_Circ,Cin.Propagation_mat_array(pp))*Phase_adjust;
                Field_Circ = Reflect_Mirror(Field_Circ,Cin.I_array(1));
            end
            
        end
        
    elseif Cin.type == 'folded'
        for pp = 1:Cin.Nb_mirror-1 % do one way
            Field_Circ = Propagate_E(Field_Circ,Cin.Propagation_mat_array(pp));
            Field_Circ = Reflect_Mirror(Field_Circ,Cin.I_array(pp+1));
        end
        for pp=Cin.Nb_mirror-1:-1:1 % and do the round trip
            Field_Circ = Propagate_E(Field_Circ,Cin.Propagation_mat_array(pp));
            if pp == 1 % last propagation
                Field_Circ = Field_Circ*Phase_adjust;
            end
            Field_Circ = Reflect_Mirror(Field_Circ,Cin.I_array(pp));
        end
    else
        error('Cavity_resonance_phase(): wrong cavity type!')
    end
    
    Phase_adjust = Phase_adjust * exp(-1i*angle(Calculate_Overlap(Field_Circ,Field_total)));
    %E_plot(Field_total); pause(0.2)
    %Calculate_power(Field_total)
end




% Then find the round trip to make the eigen mode on resonnance

Field_before = Field_total;
Field_Circ = Field_total;

if Cin.type == 'ring'
    for pp=1:Cin.Nb_mirror
        
        if pp ~= Cin.Nb_mirror % check we are not at the last iteration
            Field_Circ = Propagate_E(Field_Circ,Cin.Propagation_mat_array(pp));
            Field_Circ = Reflect_Mirror(Field_Circ,Cin.I_array(pp+1));
        else % we are at the last iteration
            Field_Circ = Propagate_E(Field_Circ,Cin.Propagation_mat_array(pp));
            Field_Circ = Reflect_Mirror(Field_Circ,Cin.I_array(1));
        end
        
    end
    
elseif Cin.type == 'folded'
    for pp = 1:Cin.Nb_mirror-1 % do one way
        Field_Circ = Propagate_E(Field_Circ,Cin.Propagation_mat_array(pp));
        Field_Circ = Reflect_Mirror(Field_Circ,Cin.I_array(pp+1));
    end
    for pp=Cin.Nb_mirror-1:-1:1 % and do the round trip
        Field_Circ = Propagate_E(Field_Circ,Cin.Propagation_mat_array(pp));
        Field_Circ = Reflect_Mirror(Field_Circ,Cin.I_array(pp));
    end
else
    error('Cavity_resonance_phase(): wrong cavity type!')
end

% Field_before.Fit_TEM00;
% Field_Circ.Fit_TEM00;
% angle(Calculate_Overlap(Field_Circ,Field_before))
Cout.Resonance_phase = exp(-1i* angle(Calculate_Overlap(Field_Circ,Field_before)));

%----------------------------------------------------------------------------------------
% Found an approximation for the resonanting field in the cavity

RT_loss = 1;
for pp=1:Cin.Nb_mirror
    RT_loss = RT_loss * Cin.I_array(pp).r;
end

Cavity_gain = (1-Cin.I_array(1).r^2)/ (1 - RT_loss)^2;
% do the overlap
Coeff_over = Calculate_Overlap(Field_Circ,Field_in2);
Pcirc = abs(Coeff_over).^2 * Cavity_gain * Calculate_Power(Cin.Laser_in);
Cout.Field_reso_guess = Normalise_E(Field_Circ,'Power',Pcirc);

disp(['Found the phase for resonance in cavity ' inputname(1)])

end