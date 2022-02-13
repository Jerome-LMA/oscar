function [] = propagate_field(obj)
%  propagate_field() Propagates and stores the transient field in
%  a cavity.
% This procedure take for argument an instance of the class Cavity1 and
% return the same object with all the transient field stored in a matrix.
% That is used to later scan the cavity.



% Find the number of iteration for the cavity. We want an accuracy of 1%

% Round trip amplitude gain for the field (= 1 - round trip loss)
RT_loss = obj.i_input.r*obj.i_end.r;

% Have to solve RT_loss^num_iter < 0.5*accuracy
num_iter = log(0.5*0.01)/(log(RT_loss));
num_iter = round(num_iter);

if (num_iter > obj.cavity_scan_param(4))
    fprintf(' Ideal number of iteration to scan the cavity: %i \n',num_iter)
    fprintf(' This is too much, number of iteration brings down to: %i \n',obj.cavity_scan_param(4))
    num_iter = obj.cavity_scan_param(4);
end

% If we do not start the beam on the input mirror, let pass through it
% first
if ~obj.laser_start_on_input
    if isa(obj.i_input, 'Interface')
        obj.laser_in =  Change_E_n(obj.laser_in,obj.i_input.n2);
    end
    obj.laser_in = Transmit_Reflect_Optic(obj.laser_in,obj.i_input);
else
    obj.laser_in =  obj.laser_in * (sqrt(1-obj.i_input.r^2));
end


Grid_num_point = obj.laser_in.Grid.Num_point;

if isgpuarray(obj.laser_in.Field)
    Total_Field = complex(zeros(Grid_num_point,Grid_num_point,num_iter,'gpuArray'));
else
    Total_Field = complex(zeros(Grid_num_point,Grid_num_point,num_iter,'double'));
end

Field_Circ = obj.laser_in;
Total_Field(:,:,1) = Field_Circ.Field;

for q = 2:num_iter
    Field_Circ = Propagate_E(Field_Circ,obj.propagation_mat);
    Field_Circ = reflect_mirror(Field_Circ,obj.i_end);
    Field_Circ = Propagate_E(Field_Circ,obj.propagation_mat);
    Field_Circ = reflect_mirror(Field_Circ,obj.i_input);
    Total_Field(:,:,q) = Field_Circ.Field;
end


obj.cavity_scan_all_field = Total_Field;

end
