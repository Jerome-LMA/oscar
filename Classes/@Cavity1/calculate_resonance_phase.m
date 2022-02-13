function [] = calculate_resonance_phase(obj, varargin)
%  Cout = Cavity_resonance_phase(Cin) find the resonance of a cavity
% This procedure can find the resonance of the cavity by finding the
% suitable round trip phase shift to bring the circulating field on resonance.

p  = inputParser;

% Check if the first argument is an interface
p.addRequired('Cin', @(x)isa(x, 'Cavity1'));

% Display or not the end
p.addParameter('verbose',true,@(x)isa(x,'logical'));

p.parse(obj,varargin{:})

% Normalise the input to 1W
obj.laser_in = Normalise_E(obj.laser_in);

if ~obj.laser_start_on_input
    if isa(obj.i_input, 'Interface')
        obj.laser_in =  Change_E_n(obj.laser_in,obj.i_input.n2);
    end
    Field_Circ1 = Transmit_Reflect_Optic(obj.laser_in,obj.i_input);
else
    Field_Circ1 = obj.laser_in;
end


num_iter = obj.cavity_phase_param;

Field_total = Normalise_E(Field_Circ1,0);
Phase_adjust =1;
Field_Circ = Field_Circ1;

% Calculate the pseudo eigen mode in the cavity
for q = 1:num_iter
    
    Field_total = Field_total + Field_Circ;
    
    % Do a round trip
    Field_Circ = Propagate_E(Field_Circ,obj.propagation_mat);
    Field_Circ = reflect_mirror(Field_Circ,obj.i_end);
    Field_Circ = Propagate_E(Field_Circ,obj.propagation_mat)*Phase_adjust;
    Field_Circ = reflect_mirror(Field_Circ,obj.i_input);
    
    Phase_adjust = Phase_adjust * exp(-1i*angle(Calculate_Overlap(Field_Circ,Field_total)));
    %     Calculate_power(Field_total)
end

% Then find the round trip to make the eigen mode on resonnance

Field_before = Field_total;
Field_Circ = Field_total;

Field_Circ = Propagate_E(Field_Circ,obj.propagation_mat);
Field_Circ = reflect_mirror(Field_Circ,obj.i_end);
%angle(Calculate_Overlap(Field_Circ,Field_before))
Field_Circ = Propagate_E(Field_Circ,obj.propagation_mat);
Field_Circ = reflect_mirror(Field_Circ,obj.i_input);

% Field_before.Fit_TEM00;
% Field_Circ.Fit_TEM00;
% angle(Calculate_Overlap(Field_Circ,Field_before))
obj.resonance_phase = exp(-1i* angle(Calculate_Overlap(Field_Circ,Field_before)));

if p.Results.verbose
    disp(['Found the phase for resonance in cavity: ' inputname(1)])
end

% Found an approximation for the resonanting field in the cavity

Cavity_gain = (1-obj.i_input.r^2)/ (1 - obj.i_input.r * obj.i_end.r)^2;
% do the overlap
Coeff_over = Calculate_Overlap(Field_Circ, Field_Circ1);
Pcirc = abs(Coeff_over).^2 * Cavity_gain; % For 1W input power
obj.field_reso_guess = Normalise_E(Field_Circ, Pcirc);

% Suppose that the SB are in antiresonance

if obj.field_reso_guess.Nb_Pair_SB
    Cavity_gain_SB = (1-obj.i_input.r^2)/ (1 + obj.i_input.r * obj.i_end.r)^2;
    Pcirc_SB = abs(Coeff_over).^2 * Cavity_gain_SB;
    
    tmp_field_SB = Normalise_E(Field_Circ,Pcirc_SB);
    
    for ii=1:obj.field_reso_guess.Nb_Pair_SB
        obj.field_reso_guess.SB(ii).Field_lower = tmp_field_SB.Field;
        obj.field_reso_guess.SB(ii).Field_upper = tmp_field_SB.Field;
    end
end

end