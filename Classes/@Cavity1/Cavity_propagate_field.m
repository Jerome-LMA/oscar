function Cout = Cavity_propagate_field(Cin)
%  Cavity_propagate_field(Cin) Propagates and stores the transient field in
%  a cavity.
% This procedure take for argument an instance of the class Cavity1 and
% return the same object with all the transient field stored in a matrix.
% That is used to later scan the cavity.



% Find the number of iteration for the cavity. We want an accuracy of 1%

% Round trip amplitude gain for the field (= 1 - round trip loss)
RT_loss = Cin.I_input.r*Cin.I_end.r;

% Have to solve RT_loss^num_iter < 0.5*accuracy
num_iter = log(0.5*0.01)/(log(RT_loss));
num_iter = round(num_iter);

if (num_iter > Cin.Cavity_scan_param(4))
    fprintf(' Ideal number of iteration to scan the cavity; %i \n',num_iter)
    fprintf(' This is too much, number of iteration brings down to; %i \n',Cin.Cavity_scan_param(4))
    num_iter = Cin.Cavity_scan_param(4);
end

Cout = Cin;

% If we do not start the beam on the input beam, let pass through it
% first
if ~Cin.Laser_start_on_input
    if isa(Cin.I_input, 'Interface')
        Cin.Laser_in =  Change_E_n(Cin.Laser_in,Cin.I_input.n2);
    end
    Cin.Laser_in = Transmit_Reflect_Optic(Cin.Laser_in,Cin.I_input);
else
    Cin.Laser_in =  Cin.Laser_in * (sqrt(1-Cin.I_input.r^2));
end


Grid_num_point = Cin.Laser_in.Grid.Num_point;

Total_Field = complex(zeros(Grid_num_point,Grid_num_point,num_iter,'double'));

Field_Circ = Cin.Laser_in;
Total_Field(:,:,1) = Field_Circ.Field;

for q = 2:num_iter
    Field_Circ = Propagate_E(Field_Circ,Cin.Propagation_mat);
    Field_Circ = Reflect_Mirror(Field_Circ,Cin.I_end);
    Field_Circ = Propagate_E(Field_Circ,Cin.Propagation_mat);
    Field_Circ = Reflect_Mirror(Field_Circ,Cin.I_input);
    Total_Field(:,:,q) = Field_Circ.Field;
end


Cout.Cavity_scan_all_field = Total_Field;

end
