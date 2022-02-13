function Cout = Calculate_Fields(Cin)
% Cout = Calculate_fields(Cin) calculate the circulating, reflected and transmitted fields
% Function used to calculated the fields inside the cavity. The laser beam must be defined outside the cavity in order to calculate the reflected field.
% Stored in Cout.field_trans(n) the transmitted field after the nth mirror.


if isempty(Cin.resonance_phase)
    error(['Calculate_fields(' inputname(1) '): The resonance position must be calculated first'])
end

if Cin.laser_start_on_input
    error(['Calculate_fields(' inputname(1) '): To calculate the reflected beam, the beam must be defined outside the cavity, set laser_start_on_input = false'])
end

Cout = Cin;

Cout.field_trans = E_Field.empty(Cin.Nb_mirror,0);
Cout.field_trans(1) = Normalise_E(Cin.laser_in,0);% But an empty field here
Cout.field_trans(1) = Change_E_n(Cout.field_trans(1),Cin.I_array(1).n2);
% Calculate the number of iteration to reach the steady state
Accuracy = 0.0001;

RT_loss = 1;
for pp=1:Cin.Nb_mirror
    RT_loss = RT_loss * Cin.I_array(pp).r;
end
% Have to solve RT_loss^num_iter < 0.5*accuracy
num_iter = log(0.5*Accuracy)/(log(RT_loss));
num_iter = round(num_iter);

Field_total = Normalise_E(Cin.laser_in,0);
Buildup_power = zeros(1,num_iter);

% The laser starts outside the input mirror, change n from 1 to mirror
% substrate refractive index
Field_in =  Change_E_n(Cin.laser_in,Cin.I_array(1).n2);
[field_transient,field_reflec] = Transmit_Reflect_Interface(Field_in,Cin.I_array(1));


for q = 1:num_iter
    
    Field_total = Field_total + field_transient;
    Buildup_power(q) = calculate_power(Field_total);
    
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
    
end

Cout.field_circ = Field_total;

%------------------------------------------------------------------
% Calculate the transmitted and reflected field

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
        end
        
    end
    
elseif Cin.type == 'folded'
    
    for pp = 1:Cin.Nb_mirror-1 % do one way
        field_transient = Propagate_E(field_transient,Cin.propagation_mat_array(pp));
        [tmp_trans,field_transient] = Transmit_Reflect_Interface(field_transient,Cin.I_array(pp+1));
        Cout.field_trans(pp+1) = tmp_trans;
        Cout.field_trans(pp+1) = Change_E_n(Cout.field_trans(pp+1),Cin.I_array(pp+1).n1);
    end
    for pp=Cin.Nb_mirror-1:-1:1 % and do the round trip
        field_transient = Propagate_E(field_transient,Cin.propagation_mat_array(pp));
        if pp == 1 % last propagation
            field_transient = field_transient*Cin.resonance_phase;
        end
        [tmp_trans,field_transient] = Transmit_Reflect_Interface(field_transient,Cin.I_array(pp));
        if pp ~= 1 % last propagation
            tmp_trans = Change_E_n(tmp_trans,Cin.I_array(pp).n1);
        end       
        Cout.field_trans(pp) = Cout.field_trans(pp) + tmp_trans;
    end
end

Cout.field_ref = field_reflec + Cout.field_trans(1);
Cout.field_ref = Change_E_n(Cout.field_ref,Cin.I_array(1).n1);



%-------------------------------------------------------------------

end