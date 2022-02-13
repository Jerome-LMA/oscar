function [Field_total,RTL] = get_info_ac(Cin)

% get_info(C1) calculate some informations about the cavity
% get_info(C1) calculate fields in the cavity, power buildup, diffraction loss
% Do not calculate the reflected field use the function 'Calculate_fields' instead.

if isempty(Cin.resonance_phase)
    error(['get_info(' inputname(1) '): The resonance position must be calculated first'])
end

if ~Cin.laser_start_on_input
    if isa(Cin.i_input, 'Interface')
        Cin.laser_in =  Change_E_n(Cin.laser_in,Cin.i_input.n2);
    end
    
    Cin.laser_in = Transmit_Reflect_Optic(Cin.laser_in,Cin.i_input,'AR');
else
    Cin.laser_in =  Cin.laser_in * Cin.i_input.t;
end

% Define the overlap function:
Raw_overlap = @(x,y) (sum(conj(x).*y, 'all') );


% Calculate the number of iteration to reach the steady state for a given
% accurary
Accuracy = 1E-12;

Field_in = Cin.laser_in;

ii = 0;
error_P = 1;
E1 = Cin.field_reso_guess;

% Find first D1 = E1 - A E1
E1_circ = Propagate_E(E1,Cin.propagation_mat);
E1_circ = reflect_mirror(E1_circ,Cin.i_end);
E1_circ = Propagate_E(E1_circ,Cin.propagation_mat);

E1_circ = E1_circ * Cin.resonance_phase;
E1_circ = reflect_mirror(E1_circ,Cin.i_input);

D1 = E1 -  E1_circ;

% On now we have E1, D1 that all we need

while error_P > Accuracy;
    ii = ii +1;
    E_SR_2 = E1 - D1 + Field_in;
    
    % Do a round trip for E_SR_2
    E_SR_2_circ = Propagate_E(E_SR_2,Cin.propagation_mat);
    E_SR_2_circ = reflect_mirror(E_SR_2_circ,Cin.i_end);
    E_SR_2_circ = Propagate_E(E_SR_2_circ,Cin.propagation_mat);
    
    E_SR_2_circ = E_SR_2_circ * Cin.resonance_phase;
    E_SR_2_circ = reflect_mirror(E_SR_2_circ,Cin.i_input);
    
    D_SR_2 = E_SR_2 - E_SR_2_circ;
    
    % Find the best coefficient a and b
    
    M(1,1) = Raw_overlap(D1.Field,D1.Field);
    M(1,2) = Raw_overlap(D1.Field,D_SR_2.Field);
    M(2,1) = conj( M(1,2) );
    M(2,2) = Raw_overlap(D_SR_2.Field,D_SR_2.Field);
    
    A(1,1) = Raw_overlap(D1.Field,Field_in.Field);
    A(2,1) = Raw_overlap(D_SR_2.Field,Field_in.Field);
    
    c = M\A;
    
    E2 = c(1)*E1 + c(2)*E_SR_2;
    
    % Calculate D2 now
    D2 = E2 - ( c(1)*(E1 - D1) + c(2)*E_SR_2_circ );
    
    error_P = Calculate_power(E2 - E1) / Calculate_power(E1);
    
    E1 = E2;
    D1 = D2;
    
end

Field_total = E1;

% Calculate the diffraction loss
field_tmp = Field_total;
field_tmp = Normalise_E(field_tmp);
field_tmp = Propagate_E(field_tmp,Cin.propagation_mat);
field_tmp = reflect_mirror(field_tmp,Cin.i_end,'Ref',1);
field_tmp = Propagate_E(field_tmp,Cin.propagation_mat);
field_tmp = reflect_mirror(field_tmp,Cin.i_input,'Ref',1);

Cavity_loss = (1 - Calculate_power(field_tmp));
RTL = Cavity_loss*1E6;

disp([' ---- Display results for cavity ' inputname(1)   ' -----'])
fprintf(' Round trip diffraction loss: %g [ppm] \n',Cavity_loss*1E6)

fprintf(' Circulating power: %g [W] \n',Calculate_power(Field_total))

% Find if the mode is a TEM00 or not
[~, m n] = Read_mode_name(Cin.laser_in.Mode_name);

if (m==0) && (n==0)
    [Beam_radius, ~] = Fit_TEM00(Propagate_E(Field_total,Cin.Length));
    fprintf(' Size of the beam on the end mirror: %g [m] \n',Beam_radius)
    
    [Beam_radius, Beam_RofC] = Fit_TEM00(Field_total);
    fprintf(' Size of the beam on the input mirror: %g [m] \n',Beam_radius)
    
    
    % Calculate the q parameter of the beam on the input mirror
    Cavity_q = 1/(1/Beam_RofC - 1i*Cin.laser_in.Wavelength/(pi*Beam_radius^2));
    Cavity_waist_size = sqrt(( (imag(Cavity_q))*Cin.laser_in.Wavelength/pi ));
    Arm.waist_position = real(Cavity_q);
    
    fprintf(' Size of the cavity waist: %g [m] \n',Cavity_waist_size)
    fprintf(' Distance of the cavity waist from the input mirror: %g [m] \n',Arm.waist_position)
end

end