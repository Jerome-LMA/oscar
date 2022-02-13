function Cout = calculate_fields_ac(Cin,varargin)
% Cout = calculate_fields_ac(Cin) calculate the circulating, reflected and transmitted fields
% Function used to calculated the fields inside the cavity. The laser beam must be defined outside the cavity in order to calculate the reflected field.
% Use the accelerated convergence scheme

p  = inputParser;
p.FunctionName = 'Calculate fields inside the cavity';

% Check if the first argument is an interface
p.addRequired('Cin', @(x)isa(x, 'CavityN'));

% Check if the resolution of the grid if given
p.addParameter('accuracy',[],@(x)isnumeric(x) && x>0);

% Check if the resolution of the grid if given
p.addParameter('iter',[],@(x)isnumeric(x) && x>0);

p.parse(Cin,varargin{:})

if isempty(Cin.resonance_phase)
    error(['Calculate_fields(' inputname(1) '): The resonance position must be calculated first'])
end

if Cin.laser_start_on_input
    error(['Calculate_fields(' inputname(1) '): To calculate the reflected beam, the beam must be defined outside the cavity, set laser_start_on_input = false'])
end

Cout = Cin;

if ~isempty(p.Results.accuracy)
    Accuracy = p.Results.accuracy;
else
    Accuracy = 1E-12;
end

% Keep some space
Cout.field_trans = E_Field.empty(Cin.Nb_mirror,0);
Cout.field_trans(1) = Normalise_E(Cin.laser_in,0);% But an empty field here

% Define the overlap function:
Raw_overlap = @(x,y) (sum(conj(x).*y, 'all') );

% The laser starts outside the input mirror, change n from 1 to mirror
% substrate refractive index

Field_in =  Change_E_n(Cin.laser_in,Cin.I_array(1).n2);
[Field_in , field_reflec] = Transmit_Reflect_Interface(Field_in,Cin.I_array(1));

%ii = 0;
error_P = 1;
E1 = Cin.field_reso_guess * sqrt(calculate_power(Cin.laser_in));  % The 'guess' field was calculated for 1W of input power, so it has to be normalised according to the current incident power

if ~E1.Nb_Pair_SB % if there is no SB
    
    % Find first D1 = E1 - A E1
    field_transient = E1;
    
    if Cin.type == 'ring'
        for pp=1:Cin.Nb_mirror
            
            if pp ~= Cin.Nb_mirror % check we are not at the last iteration
                field_transient = Propagate_E(field_transient,Cin.propagation_mat_array(pp));
                field_transient = reflect_mirror(field_transient,Cin.I_array(pp+1));
            else % we are at the last iteration
                field_transient = Propagate_E(field_transient,Cin.propagation_mat_array(pp))* Cin.resonance_phase;
                field_transient = reflect_mirror(field_transient,Cin.I_array(1));
            end
            
        end
        
    elseif Cin.type == 'folded'
        for pp = 1:Cin.Nb_mirror-1 % do one way
            field_transient = Propagate_E(field_transient,Cin.propagation_mat_array(pp));
            field_transient = reflect_mirror(field_transient,Cin.I_array(pp+1));
        end
        for pp=Cin.Nb_mirror-1:-1:1 % and do the round trip
            field_transient = Propagate_E(field_transient,Cin.propagation_mat_array(pp));
            if pp == 1 % last propagation
                field_transient = field_transient*Cin.resonance_phase;
            end
            field_transient = reflect_mirror(field_transient,Cin.I_array(pp));
        end
    end
    
    
    D1 = E1 -  field_transient;
    
    % On now we have E1, D1 that all we need
    
    while error_P > Accuracy
        %  ii = ii +1
        E_SR_2 = E1 - D1 + Field_in;
        
        % Do a round trip for E_SR_2
        field_transient = E_SR_2;
        
        if Cin.type == 'ring'
            for pp=1:Cin.Nb_mirror
                
                if pp ~= Cin.Nb_mirror % check we are not at the last iteration
                    field_transient = Propagate_E(field_transient,Cin.propagation_mat_array(pp));
                    field_transient = reflect_mirror(field_transient,Cin.I_array(pp+1));
                else % we are at the last iteration
                    field_transient = Propagate_E(field_transient,Cin.propagation_mat_array(pp))* Cin.resonance_phase;
                    field_transient = reflect_mirror(field_transient,Cin.I_array(1));
                end
                
            end
            
        elseif Cin.type == 'folded'
            for pp = 1:Cin.Nb_mirror-1 % do one way
                field_transient = Propagate_E(field_transient,Cin.propagation_mat_array(pp));
                field_transient = reflect_mirror(field_transient,Cin.I_array(pp+1));
            end
            for pp=Cin.Nb_mirror-1:-1:1 % and do the round trip
                field_transient = Propagate_E(field_transient,Cin.propagation_mat_array(pp));
                if pp == 1 % last propagation
                    field_transient = field_transient*Cin.resonance_phase;
                end
                field_transient = reflect_mirror(field_transient,Cin.I_array(pp));
            end
        end
        
        E_SR_2_circ = field_transient;
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
        
        error_P = calculate_power(E2 - E1) / calculate_power(E1);
        
        E1 = E2;
        D1 = D2;
        
    end
    
    Field_total = E1;
    Cout.field_circ = Field_total; % circulating field found
    
else
    error('calculate_fields_ac(): Calcuations for the sidebands not yet implemented: contact the developers if it is needed')
end

%---------------------------------------------------------------------------------------------------------------------------------
% Ok we got now the circulating field
field_transient = Field_total;

if Cin.type == 'ring'
    for pp=1:Cin.Nb_mirror
        
        if pp ~= Cin.Nb_mirror % check we are not at the last iteration
            field_transient = Propagate_E(field_transient,Cin.propagation_mat_array(pp));
            [Cout.field_trans(pp+1),field_transient] = Transmit_Reflect_Interface(field_transient,Cin.I_array(pp+1));
            Cout.field_trans(pp+1) = Change_E_n(Cout.field_trans(pp+1),Cin.I_array(pp+1).n1);
        else % we are at the last iteration
            field_transient = Propagate_E(field_transient,Cin.propagation_mat_array(pp))* Cin.resonance_phase;
            [Cout.field_trans(1),field_transient] = Transmit_Reflect_Interface(field_transient,Cin.I_array(1));
            Cout.field_trans(1) = Change_E_n(Cout.field_trans(1),Cin.I_array(1).n1);
        end
    end
    
elseif Cin.type == 'folded'
    
    for pp = 1:Cin.Nb_mirror-1 % do one way
        %calculate_power(field_transient)
        field_transient = Propagate_E(field_transient,Cin.propagation_mat_array(pp));
        [tmp_trans , field_transient] = Transmit_Reflect_Interface(field_transient,Cin.I_array(pp+1));
        Cout.field_trans(pp+1) = tmp_trans;
        Cout.field_trans(pp+1) = Change_E_n(Cout.field_trans(pp+1),Cin.I_array(pp+1).n1);
    end
    
    % Add the round trip phase adjustement
    field_transient = field_transient*Cin.resonance_phase;
    
    for pp=Cin.Nb_mirror-1:-1:1 % and do the round trip
        %calculate_power(field_transient)
        field_transient = Propagate_E(field_transient,Cin.propagation_mat_array(pp));
        [tmp_trans,field_transient] = Transmit_Reflect_Interface(field_transient,Cin.I_array(pp));
        tmp_trans = Change_E_n(tmp_trans,Cin.I_array(pp).n1);
        Cout.field_trans(pp) = Cout.field_trans(pp) + tmp_trans;
    end
end

field_reflec = Change_E_n(field_reflec,Cin.I_array(1).n1);
Cout.field_ref = field_reflec + Cout.field_trans(1);

%-------------------------------------------------------------------

end