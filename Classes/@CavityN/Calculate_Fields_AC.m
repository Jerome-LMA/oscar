function Cout = Calculate_Fields_AC(Cin,varargin)
% Cout = Calculate_fields_AC(Cin) calculate the circulating, reflected and transmitted fields
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

if isempty(Cin.Resonance_phase)
    error(['Calculate_fields(' inputname(1) '): The resonance position must be calculated first'])
end

if Cin.Laser_start_on_input
    error(['Calculate_fields(' inputname(1) '): To calculate the reflected beam, the beam must be defined outside the cavity, set Laser_start_on_input = false'])
end

Cout = Cin;

if ~isempty(p.Results.accuracy)
    Accuracy = p.Results.accuracy;
else
    Accuracy = 1E-12;
end

% Keep some space
Cout.Field_trans = E_Field.empty(Cin.Nb_mirror,0);
Cout.Field_trans(1) = Normalise_E(Cin.Laser_in,0);% But an empty field here

% Define the overlap function:
Raw_overlap = @(x,y) (sum(conj(x).*y, 'all') );

% The laser starts outside the input mirror, change n from 1 to mirror
% substrate refractive index

Field_in =  Change_E_n(Cin.Laser_in,Cin.I_array(1).n2);
[Field_in , Field_reflec] = Transmit_Reflect_Interface(Field_in,Cin.I_array(1));

%ii = 0;
error_P = 1;
E1 = Cin.Field_reso_guess * sqrt(Calculate_Power(Cin.Laser_in));  % The 'guess' field was calculated for 1W of input power, so it has to be normalised according to the current incident power

if ~E1.Nb_Pair_SB % if there is no SB
    
    % Find first D1 = E1 - A E1
    Field_transient = E1;
    
    if Cin.type == 'ring'
        for pp=1:Cin.Nb_mirror
            
            if pp ~= Cin.Nb_mirror % check we are not at the last iteration
                Field_transient = Propagate_E(Field_transient,Cin.Propagation_mat_array(pp));
                Field_transient = Reflect_Mirror(Field_transient,Cin.I_array(pp+1));
            else % we are at the last iteration
                Field_transient = Propagate_E(Field_transient,Cin.Propagation_mat_array(pp))* Cin.Resonance_phase;
                Field_transient = Reflect_Mirror(Field_transient,Cin.I_array(1));
            end
            
        end
        
    elseif Cin.type == 'folded'
        for pp = 1:Cin.Nb_mirror-1 % do one way
            Field_transient = Propagate_E(Field_transient,Cin.Propagation_mat_array(pp));
            Field_transient = Reflect_Mirror(Field_transient,Cin.I_array(pp+1));
        end
        for pp=Cin.Nb_mirror-1:-1:1 % and do the round trip
            Field_transient = Propagate_E(Field_transient,Cin.Propagation_mat_array(pp));
            if pp == 1 % last propagation
                Field_transient = Field_transient*Cin.Resonance_phase;
            end
            Field_transient = Reflect_Mirror(Field_transient,Cin.I_array(pp));
        end
    end
    
    
    D1 = E1 -  Field_transient;
    
    % On now we have E1, D1 that all we need
    
    while error_P > Accuracy
        %  ii = ii +1
        E_SR_2 = E1 - D1 + Field_in;
        
        % Do a round trip for E_SR_2
        Field_transient = E_SR_2;
        
        if Cin.type == 'ring'
            for pp=1:Cin.Nb_mirror
                
                if pp ~= Cin.Nb_mirror % check we are not at the last iteration
                    Field_transient = Propagate_E(Field_transient,Cin.Propagation_mat_array(pp));
                    Field_transient = Reflect_Mirror(Field_transient,Cin.I_array(pp+1));
                else % we are at the last iteration
                    Field_transient = Propagate_E(Field_transient,Cin.Propagation_mat_array(pp))* Cin.Resonance_phase;
                    Field_transient = Reflect_Mirror(Field_transient,Cin.I_array(1));
                end
                
            end
            
        elseif Cin.type == 'folded'
            for pp = 1:Cin.Nb_mirror-1 % do one way
                Field_transient = Propagate_E(Field_transient,Cin.Propagation_mat_array(pp));
                Field_transient = Reflect_Mirror(Field_transient,Cin.I_array(pp+1));
            end
            for pp=Cin.Nb_mirror-1:-1:1 % and do the round trip
                Field_transient = Propagate_E(Field_transient,Cin.Propagation_mat_array(pp));
                if pp == 1 % last propagation
                    Field_transient = Field_transient*Cin.Resonance_phase;
                end
                Field_transient = Reflect_Mirror(Field_transient,Cin.I_array(pp));
            end
        end
        
        E_SR_2_circ = Field_transient;
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
        
        error_P = Calculate_Power(E2 - E1) / Calculate_Power(E1);
        
        E1 = E2;
        D1 = D2;
        
    end
    
    Field_total = E1;
    Cout.Field_circ = Field_total; % circulating field found
    
else
    error('Calculate_fields_AC(): Calcuations for the sidebands not yet implemented: contact the developers if it is needed')
end

%---------------------------------------------------------------------------------------------------------------------------------
% Ok we got now the circulating field
Field_transient = Field_total;

if Cin.type == 'ring'
    for pp=1:Cin.Nb_mirror
        
        if pp ~= Cin.Nb_mirror % check we are not at the last iteration
            Field_transient = Propagate_E(Field_transient,Cin.Propagation_mat_array(pp));
            [Cout.Field_trans(pp+1),Field_transient] = Transmit_Reflect_Interface(Field_transient,Cin.I_array(pp+1));
            Cout.Field_trans(pp+1) = Change_E_n(Cout.Field_trans(pp+1),Cin.I_array(pp+1).n1);
        else % we are at the last iteration
            Field_transient = Propagate_E(Field_transient,Cin.Propagation_mat_array(pp))* Cin.Resonance_phase;
            [Cout.Field_trans(1),Field_transient] = Transmit_Reflect_Interface(Field_transient,Cin.I_array(1));
            Cout.Field_trans(1) = Change_E_n(Cout.Field_trans(1),Cin.I_array(1).n1);
        end
    end
    
elseif Cin.type == 'folded'
    
    for pp = 1:Cin.Nb_mirror-1 % do one way
        %Calculate_Power(Field_transient)
        Field_transient = Propagate_E(Field_transient,Cin.Propagation_mat_array(pp));
        [tmp_trans , Field_transient] = Transmit_Reflect_Interface(Field_transient,Cin.I_array(pp+1));
        Cout.Field_trans(pp+1) = tmp_trans;
        Cout.Field_trans(pp+1) = Change_E_n(Cout.Field_trans(pp+1),Cin.I_array(pp+1).n1);
    end
    
    % Add the round trip phase adjustement
    Field_transient = Field_transient*Cin.Resonance_phase;
    
    for pp=Cin.Nb_mirror-1:-1:1 % and do the round trip
        %Calculate_Power(Field_transient)
        Field_transient = Propagate_E(Field_transient,Cin.Propagation_mat_array(pp));
        [tmp_trans,Field_transient] = Transmit_Reflect_Interface(Field_transient,Cin.I_array(pp));
        tmp_trans = Change_E_n(tmp_trans,Cin.I_array(pp).n1);
        Cout.Field_trans(pp) = Cout.Field_trans(pp) + tmp_trans;
    end
end

Field_reflec = Change_E_n(Field_reflec,Cin.I_array(1).n1);
Cout.Field_ref = Field_reflec + Cout.Field_trans(1);

%-------------------------------------------------------------------

end