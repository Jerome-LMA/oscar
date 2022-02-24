function Cout = Cavity_propagate_field(Cin)
%  Cavity_propagate_field(Cin) Propagates and stores the transient field in
%  a cavity.
% This procedure take for argument an instance of the class Cavity1 and
% return the same object with all the transient field stored in a matrix.
% That is used to later scan the cavity.



% Find the number of iteration for the cavity. We want an accuracy of 1%

% Round trip amplitude gain for the field (= 1 - round trip loss)
RT_loss = 1;
for pp=1:Cin.Nb_mirror
    RT_loss = RT_loss * Cin.I_array(pp).r;
end

% Have to solve RT_loss^num_iter < 0.5*accuracy
num_iter = log(0.5*0.01)/(log(RT_loss));
num_iter = round(num_iter);

if (num_iter > 200)
    fprintf(' Number of iteration to scan the cavity; %i \n',num_iter)
end

Cout = Cin;

% If we do not start the beam on the input beam, let pass through it
% first
if ~Cin.Laser_start_on_input
    Field_in =  Change_E_n(Cin.Laser_in,Cin.I_array(1).n2);
    Field_in = Transmit_Reflect_Interface(Field_in,Cin.I_array(1));
    Field_Circ = Field_in;
else
    Field_Circ = Cin.Laser_in * Cin.I_array(1).t;
end

Grid_num_point = Cin.Laser_in.Grid.Num_point;
Total_Field = complex(zeros(Grid_num_point,Grid_num_point,num_iter,'double'));

Total_Field(:,:,1) = Field_Circ.Field;

for q = 2:num_iter
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
    end
    
    Total_Field(:,:,q) = Field_Circ.Field;
end



Cout.Cavity_scan_all_field = Total_Field;

end
