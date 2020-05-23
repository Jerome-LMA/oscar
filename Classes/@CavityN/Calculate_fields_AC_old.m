function Cout = Calculate_fields_AC(Cin,varargin)
% Cout = Calculate_fields_AC(Cin) calculate the circulating, reflected and transmitted fields
% Function used to calculated the fields inside the cavity. The laser beam must be defined outside the cavity in order to calculate the reflected field.
% Stored in Cout.Field_trans(n) the transmitted field after the nth mirror.
% Use the accelerated convergence scheme

p  = inputParser;
p.FunctionName = 'Calculate fields inside the cavity';

% Check if the first argument is an interface
p.addRequired('Cin', @(x)isa(x, 'Cavity1'));

% Check if the resolution of the grid if given
p.addParamValue('accuracy',[],@(x)isnumeric(x) && x>0);

% Check if the resolution of the grid if given
p.addParamValue('iter',[],@(x)isnumeric(x) && x>0);

p.parse(Cin,varargin{:})

if isempty(Cin.Resonance_phase)
    error(['Calculate_fields(' inputname(1) '): The resonance position must be calculated first'])
end

if Cin.Laser_start_on_input
    error(['Calculate_fields(' inputname(1) '): To calculate the reflected beam, the beam must be defined outside the cavity, set Laser_start_on_input = false'])
end

% Accelerated convergence only works for the carrier. Remove the SB fields
Cin.Laser_in.Field_SBl=[];
Cin.Laser_in.Field_SBu=[];
Cin.Field_reso_guess.Field_SBl=[];
Cin.Field_reso_guess.Field_SBu=[];

Cout = Cin;

if ~isempty(p.Results.accuracy)
    Accuracy = p.Results.accuracy;
else
    Accuracy = 1E-12;
end

% Define the overlap function:
Raw_overlap = @(x,y) (sum(conj(x).*y, 'all') );

% Accelerated convergence only works for the carrier. Remove the SB fields
Cin.Laser_in.Field_SBl=[];
Cin.Laser_in.Field_SBu=[];

Cout.Field_trans = E_Field.empty(Cin.Nb_mirror,0);
Cout.Field_trans(1) = Cin.Laser_in; % Need to initiate the first value with a dummy thing, will be overwritten later.

Field_total = Normalise_E(Cin.Laser_in,0);

% The laser starts outside the input mirror, change n from 1 to mirror
% substrate refractive index
Field_in =  Change_E_n(Cin.Laser_in,Cin.I_array(1).n2);
[Field_in,Field_reflec] = Transmit_Reflect_Interface(Field_in,Cin.I_array(1));

error_P = 1;
E1 = Cin.Field_reso_guess;

% Find first D1 = E1 - A E1
E1_circ = Propagate_E(E1,Cin.Propagation_mat);
E1_circ = Reflect_mirror(E1_circ,Cin.I_end);
E1_circ = Propagate_E(E1_circ,Cin.Propagation_mat);

E1_circ = E1_circ * Cin.Resonance_phase;
E1_circ = Reflect_mirror(E1_circ,Cin.I_input);

D1 = E1 -  E1_circ;

E1_circ = E1;

    for pp=1:Cin.Nb_mirror
        if pp ~= Cin.Nb_mirror % check we are not at the last iteration
            E1_circ = Propagate_E(E1_circ,Cin.Propagation_mat_array(pp));
            E1_circ = Reflect_mirror(E1_circ,Cin.I_array(pp+1));
        else
            E1_circ = Propagate_E(E1_circ,Cin.Propagation_mat_array(pp)) * Cin.Resonance_phase;
            E1_circ = Reflect_mirror(E1_circ,Cin.I_array(1));
        end
    end

D1 = E1 -  E1_circ;

% On now we have E1, D1 that all we need

while error_P > Accuracy;
    ii = ii +1;
    E_SR_2 = E1 - D1 + Field_in;
    
    % Do a round trip for E_SR_2
    E_SR_2_circ = E_SR_2;
    
        for pp=1:Cin.Nb_mirror
        if pp ~= Cin.Nb_mirror % check we are not at the last iteration
            E_SR_2_circ = Propagate_E(E_SR_2_circ,Cin.Propagation_mat_array(pp));
            E_SR_2_circ = Reflect_mirror(E_SR_2_circ,Cin.I_array(pp+1));
        else
            E_SR_2_circ = Propagate_E(E_SR_2_circ,Cin.Propagation_mat_array(pp)) * Cin.Resonance_phase;
            E_SR_2_circ = Reflect_mirror(E_SR_2_circ,Cin.I_array(1));
        end
    end
        
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

Cout.Field_circ = E1;

%------------------------------------------------------------------
% Calculate the transmitted and reflected field

Field_transient = E1;

for pp=1:Cin.Nb_mirror
    if pp ~= Cin.Nb_mirror % check we are not at the last iteration
        Field_transient = Propagate_E(Field_transient,Cin.Propagation_mat_array(pp));
        [Cout.Field_trans(pp+1)  Field_transient] = Transmit_Reflect_Interface(Field_transient,Cin.I_array(pp+1));
    else
        Field_transient = Propagate_E(Field_transient,Cin.Propagation_mat_array(pp)) * Cin.Resonance_phase;
        [Cout.Field_trans(1)  Field_transient] = Transmit_Reflect_Interface(Field_transient,Cin.I_array(1));
    end
end

Cout.Field_ref = Field_reflec + Cout.Field_trans(1);
Cout.Field_ref =  Change_E_n(Cout.Field_ref,Cin.I_array(1).n1);
%-------------------------------------------------------------------

end