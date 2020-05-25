function Cout = Calculate_fields_AC(Cin,varargin)
% Cout = Calculate_fields_AC(Cin) calculate the circulating, reflected and transmitted fields
% Function used to calculated the fields inside the cavity. The laser beam must be defined outside the cavity in order to calculate the reflected field.
% Use the accelerated convergence scheme

p  = inputParser;
p.FunctionName = 'Calculate fields inside the cavity';

% Check if the first argument is an interface
p.addRequired('Cin', @(x)isa(x, 'Cavity1'));

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

Display_debug = false;

% Define the overlap function:
Raw_overlap = @(x,y) (sum(sum(conj(x).*y) ) );

% The laser starts outside the input mirror, change n from 1 to mirror
% substrate refractive index

if isa(Cin.I_input, 'Interface')
    Field_in =  Change_E_n(Cin.Laser_in,Cin.I_input.n2);
else
    Field_in =  Cin.Laser_in;
end

[Field_in,Field_reflec] = Transmit_Reflect_Optic(Field_in,Cin.I_input,'AR');

error_P = 1;

if isempty(Cin.Field_reso_guess) % If the cavity resonance phase has not been calculate beforehand, so no guess for the resonant field
    E1 = Field_in;
else
    E1 = Cin.Field_reso_guess * sqrt(Calculate_power(Cin.Laser_in)); % The 'guess' field was calculated for 1W of input power, so it has to be normalised according to the current incident power
end

if ~E1.Nb_Pair_SB % if there is no SB
    
    % Find first D1 = E1 - A E1
    E1_circ = Propagate_E(E1,Cin.Propagation_mat);
    E1_circ = Reflect_mirror(E1_circ,Cin.I_end);
    E1_circ = Propagate_E(E1_circ,Cin.Propagation_mat);
    
    E1_circ = E1_circ * Cin.Resonance_phase;
    E1_circ = Reflect_mirror(E1_circ,Cin.I_input);
    
    D1 = E1 -  E1_circ;
    
    % On now we have E1, D1 that all we need
    
    count_iter = 0;
    count_iter_LSB = 0;
    
    while error_P > Accuracy
        %  ii = ii +1
        E_SR_2 = E1 - D1 + Field_in;
        
        % Do a round trip for E_SR_2
        E_SR_2_circ = Propagate_E(E_SR_2,Cin.Propagation_mat);
        E_SR_2_circ = Reflect_mirror(E_SR_2_circ,Cin.I_end);
        E_SR_2_circ = Propagate_E(E_SR_2_circ,Cin.Propagation_mat);
        
        E_SR_2_circ = E_SR_2_circ * Cin.Resonance_phase;
        E_SR_2_circ = Reflect_mirror(E_SR_2_circ,Cin.I_input);
        
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
        
        count_iter = count_iter + 1;
    end
    
    if Display_debug
        fprintf('Number of iteration carrier: %i \n',count_iter)
    end
    
    Field_total = E1;
    Cout.Field_circ = Field_total; % circulating field found
    
else
    
    % SB are present, will do independently sevral iterations, one for the carrier and one for each sidebands.
    
    %---------------------------------------------------------------------------------------------------------------------------------
    % Start with the carrier
    %E1 = Cin.Field_reso_guess * sqrt(Calculate_power(Cin.Laser_in));
    
    % Find first D1 = E1 - A E1
    E1_circ = Propagate_E(E1,Cin.Propagation_mat);
    E1_circ = Reflect_mirror(E1_circ,Cin.I_end);
    E1_circ = Propagate_E(E1_circ,Cin.Propagation_mat);
    
    E1_circ = E1_circ * Cin.Resonance_phase;
    E1_circ = Reflect_mirror(E1_circ,Cin.I_input);
    
    D1 = E1 -  E1_circ;
    % On now we have E1, D1 that all we need
    
    % Initiate the figures of merit (must tend toward 0)
    error_P_carr = 1;
    error_P_LSB(1:E1.Nb_Pair_SB) = 1;
    error_P_USB(1:E1.Nb_Pair_SB) = 1;
    
    % Initiate the counters for debugging
    count_iter_car = 0;
    count_iter_LSB(1:E1.Nb_Pair_SB) = 0;
    count_iter_USB(1:E1.Nb_Pair_SB) = 0;
    
    while error_P > Accuracy
        %  ii = ii +1
        E_SR_2 = E1 - D1 + Field_in;
        
        % Do a round trip for E_SR_2
        E_SR_2_circ = Propagate_E(E_SR_2,Cin.Propagation_mat);
        E_SR_2_circ = Reflect_mirror(E_SR_2_circ,Cin.I_end);
        E_SR_2_circ = Propagate_E(E_SR_2_circ,Cin.Propagation_mat);
        
        E_SR_2_circ = E_SR_2_circ * Cin.Resonance_phase;
        E_SR_2_circ = Reflect_mirror(E_SR_2_circ,Cin.I_input);
        
        D_SR_2 = E_SR_2 - E_SR_2_circ;
        
        %-----------------------------------------------------------------------------------------------------------------
        % Find the best coefficient a and b for the carrier
        
        if error_P_carr > Accuracy
            
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
            
            error_P_carr = Calculate_power(E2 - E1) / Calculate_power(E1);
            
            E1 = E2;
            D1 = D2;
            count_iter_car = count_iter_car + 1;
        else
            error_P_carr = 0;
        end
        
        %-----------------------------------------------------------------------------------------------------------------
        % Find the best coefficient a and b for the SB, do one by one the
        % iteration for each SB fields, first lower and then upper
        
        for ii=1:E1.Nb_Pair_SB
            if error_P_LSB(ii) > Accuracy
                
                M(1,1) = Raw_overlap(D1.SB(ii).Field_lower,D1.SB(ii).Field_lower);
                M(1,2) = Raw_overlap(D1.SB(ii).Field_lower,D_SR_2.SB(ii).Field_lower);
                M(2,1) = conj( M(1,2) );
                M(2,2) = Raw_overlap(D_SR_2.SB(ii).Field_lower,D_SR_2.SB(ii).Field_lower);
                
                A(1,1) = Raw_overlap(D1.SB(ii).Field_lower,Field_in.SB(ii).Field_lower);
                A(2,1) = Raw_overlap(D_SR_2.SB(ii).Field_lower,Field_in.SB(ii).Field_lower);
                
                %c = M\A;
                c = pinv(M) * A; % To avoid error message due to singular matrix
                
                E2.SB(ii).Field_lower = c(1)*E1.SB(ii).Field_lower + c(2)*E_SR_2.SB(ii).Field_lower;
                
                % Calculate D2 now
                D2.SB(ii).Field_lower = E2.SB(ii).Field_lower - ( c(1)*(E1.SB(ii).Field_lower - D1.SB(ii).Field_lower)...
                    + c(2)*E_SR_2_circ.SB(ii).Field_lower);
                
                [pow_diff,~] = Calculate_power_SB(E2 - E1,'SB_num',ii);
                [pow_E1,~] = Calculate_power_SB(E1,'SB_num',ii);
                
                error_P_LSB(ii) = pow_diff / pow_E1;
                
                E1.SB(ii).Field_lower = E2.SB(ii).Field_lower;
                D1.SB(ii).Field_lower = D2.SB(ii).Field_lower;
                count_iter_LSB(ii) = count_iter_LSB(ii) + 1;
            else
                error_P_LSB(ii) = 0;
            end
            
            if error_P_USB(ii) > Accuracy
                
                M(1,1) = Raw_overlap(D1.SB(ii).Field_upper,D1.SB(ii).Field_upper);
                M(1,2) = Raw_overlap(D1.SB(ii).Field_upper,D_SR_2.SB(ii).Field_upper);
                M(2,1) = conj( M(1,2) );
                M(2,2) = Raw_overlap(D_SR_2.SB(ii).Field_upper,D_SR_2.SB(ii).Field_upper);
                
                A(1,1) = Raw_overlap(D1.SB(ii).Field_upper,Field_in.SB(ii).Field_upper);
                A(2,1) = Raw_overlap(D_SR_2.SB(ii).Field_upper,Field_in.SB(ii).Field_upper);
                
                %c = M\A;
                c = pinv(M) * A; % To avoid error message due to singular matrix
                
                E2.SB(ii).Field_upper = c(1)*E1.SB(ii).Field_upper + c(2)*E_SR_2.SB(ii).Field_upper;
                
                % Calculate D2 now
                D2.SB(ii).Field_upper = E2.SB(ii).Field_upper - ( c(1)*(E1.SB(ii).Field_upper - D1.SB(ii).Field_upper) + c(2)*E_SR_2_circ.SB(ii).Field_upper);
                
                [~,pow_diff] = Calculate_power_SB(E2 - E1,'SB_num',ii);
                [~,pow_E1] = Calculate_power_SB(E1,'SB_num',ii);
                
                error_P_USB(ii) = pow_diff / pow_E1;
                
                E1.SB(ii).Field_upper = E2.SB(ii).Field_upper;
                D1.SB(ii).Field_upper = D2.SB(ii).Field_upper;
                
                D1.SB(ii).Field_upper = D2.SB(ii).Field_upper;
                count_iter_USB(ii) = count_iter_USB(ii) + 1;
            else
                error_P_USB(ii) = 0;
            end      
        end
                
        error_P = max([error_P_carr error_P_LSB error_P_USB]);
    end
    
    Cout.Field_circ = E1;
    
    if Display_debug
        fprintf('Number of iteration carrier: %i \n',count_iter_car)
        fprintf('Number of iteration LSB: %i \n',count_iter_LSB)
        fprintf('Number of iteration USB: %i \n',count_iter_USB)
    end
end

%---------------------------------------------------------------------------------------------------------------------------------
% Ok we got now the 3 circulating fields (carrier + 2 sidebands)
Field_temp = Propagate_E(Cout.Field_circ,Cin.Propagation_mat);
Cout.Field_trans = Transmit_Reflect_Optic(Field_temp,Cin.I_end);

Field_temp = Propagate_E(Cout.Field_circ,Cin.Propagation_mat);
Field_temp = Reflect_mirror(Field_temp,Cin.I_end);
Field_temp = Propagate_E(Field_temp,Cin.Propagation_mat);
Field_temp = Field_temp * Cin.Resonance_phase;

Field_temp = Transmit_Reflect_Optic(Field_temp,Cin.I_input);
Cout.Field_ref = Field_reflec + Field_temp;

%
if isa(Cin.I_input, 'Interface')
    Cout.Field_ref =  Change_E_n(Cout.Field_ref,Cin.I_input.n1);
end

if isa(Cin.I_end, 'Interface')
    Cout.Field_trans =  Change_E_n(Cout.Field_trans,Cin.I_end.n1);
end

%-------------------------------------------------------------------
% Calculate the round trip loss
% First get the HR interface
if isa(Cin.I_input,'Mirror')
    Cin.I_input = Cin.I_input.I_HR;
end

if isa(Cin.I_end,'Mirror')
    Cin.I_end = Cin.I_end.I_HR;
end

field_tmp = Cout.Field_circ;
field_tmp = Normalise_E(field_tmp);
field_tmp = Propagate_E(field_tmp,Cin.Propagation_mat);
field_tmp = Reflect_mirror(field_tmp,Cin.I_end,'Ref',1);
field_tmp = Propagate_E(field_tmp,Cin.Propagation_mat);
field_tmp = Reflect_mirror(field_tmp,Cin.I_input,'Ref',1);

Cout.Loss_RTL =  (1 - Calculate_power(field_tmp));

end