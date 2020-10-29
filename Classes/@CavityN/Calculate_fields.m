function Cout = Calculate_fields(Cin)
% Cout = Calculate_fields(Cin) calculate the circulating, reflected and transmitted fields
% Function used to calculated the fields inside the cavity. The laser beam must be defined outside the cavity in order to calculate the reflected field.
% Stored in Cout.Field_trans(n) the transmitted field after the nth mirror.


if isempty(Cin.Resonance_phase)
    error(['Calculate_fields(' inputname(1) '): The resonance position must be calculated first'])
end

if Cin.Laser_start_on_input
    error(['Calculate_fields(' inputname(1) '): To calculate the reflected beam, the beam must be defined outside the cavity, set Laser_start_on_input = false'])
end

Cout = Cin;

Cout.Field_trans = E_Field.empty(Cin.Nb_mirror,0);
Cout.Field_trans(1) = Normalise_E(Cin.Laser_in,0);% But an empty field here
Cout.Field_trans(1) = Change_E_n(Cout.Field_trans(1),Cin.I_array(1).n2);
% Calculate the number of iteration to reach the steady state
Accuracy = 0.0001;

RT_loss = 1;
for pp=1:Cin.Nb_mirror
    RT_loss = RT_loss * Cin.I_array(pp).r;
end
% Have to solve RT_loss^num_iter < 0.5*accuracy
num_iter = log(0.5*Accuracy)/(log(RT_loss));
num_iter = round(num_iter);

Field_total = Normalise_E(Cin.Laser_in,0);
Buildup_power = zeros(1,num_iter);

% The laser starts outside the input mirror, change n from 1 to mirror
% substrate refractive index
Field_in =  Change_E_n(Cin.Laser_in,Cin.I_array(1).n2);
[Field_transient,Field_reflec] = Transmit_Reflect_Interface(Field_in,Cin.I_array(1));


for q = 1:num_iter
    
    Field_total = Field_total + Field_transient;
    Buildup_power(q) = Calculate_Power(Field_total);
    
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
    
end

Cout.Field_circ = Field_total;

%------------------------------------------------------------------
% Calculate the transmitted and reflected field

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
        end
        
    end
    
elseif Cin.type == 'folded'
    
    for pp = 1:Cin.Nb_mirror-1 % do one way
        Field_transient = Propagate_E(Field_transient,Cin.Propagation_mat_array(pp));
        [tmp_trans,Field_transient] = Transmit_Reflect_Interface(Field_transient,Cin.I_array(pp+1));
        Cout.Field_trans(pp+1) = tmp_trans;
        Cout.Field_trans(pp+1) = Change_E_n(Cout.Field_trans(pp+1),Cin.I_array(pp+1).n1);
    end
    for pp=Cin.Nb_mirror-1:-1:1 % and do the round trip
        Field_transient = Propagate_E(Field_transient,Cin.Propagation_mat_array(pp));
        if pp == 1 % last propagation
            Field_transient = Field_transient*Cin.Resonance_phase;
        end
        [tmp_trans,Field_transient] = Transmit_Reflect_Interface(Field_transient,Cin.I_array(pp));
        if pp ~= 1 % last propagation
            tmp_trans = Change_E_n(tmp_trans,Cin.I_array(pp).n1);
        end       
        Cout.Field_trans(pp) = Cout.Field_trans(pp) + tmp_trans;
    end
end

Cout.Field_ref = Field_reflec + Cout.Field_trans(1);
Cout.Field_ref = Change_E_n(Cout.Field_ref,Cin.I_array(1).n1);



%-------------------------------------------------------------------

end