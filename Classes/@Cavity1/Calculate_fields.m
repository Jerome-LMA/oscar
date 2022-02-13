function [] = calculate_fields(obj, varargin)
% calculate_fields: calculate the circulating, reflected and transmitted fields
% Function used to calculated the fields inside the cavity. The laser beam must be defined outside the cavity in order to calculate the reflected field.

p  = inputParser;
p.FunctionName = 'Calculate fields inside the cavity';

% Check if the resolution of the grid if given
p.addParameter('accuracy',[],@(x)isnumeric(x) && x>0);

% Check if the resolution of the grid if given
p.addParameter('iter',[],@(x)isnumeric(x) && x>0);

p.parse(obj,varargin{:})


if isempty(obj.resonance_phase)
    error(['calculate_fields(' inputname(1) '): The resonance position must be calculated first'])
end

if obj.laser_start_on_input
    error(['calculate_fields(' inputname(1) '): To calculate the reflected beam, the beam must be defined outside the cavity, set laser_start_on_input = false'])
end

if ~isempty(p.Results.accuracy)
    accuracy = p.Results.accuracy;
else
    accuracy = 0.0001;
end
% Calculate the number of iteration to reach the steady state


rt_loss = obj.i_input.r*obj.i_end.r;
% Have to solve RT_loss^num_iter < 0.5*accuracy
num_iter = log(0.5*accuracy)/(log(rt_loss));

if ~isempty(p.Results.iter) % overide the numbers of iteration
    num_iter = p.Results.iter;
end
num_iter = round(num_iter);

% The laser starts outside the input mirror, change n from 1 to mirror
% substrate refractive index

if isa(obj.i_input, 'Interface')
    Field_in =  Change_E_n(obj.laser_in,obj.i_input.n2);
else
    Field_in =  obj.laser_in;
end

[Field_transient,Field_reflec] = Transmit_Reflect_Optic(Field_in,obj.i_input,'AR');
Field_total = Normalise_E(Field_transient,0);

power_buildup = zeros(1,num_iter);

for q = 1:num_iter
    power_buildup(q) = calculate_power(Field_total);
    Field_total = Field_total + Field_transient;
    Field_transient = Propagate_E(Field_transient,obj.propagation_mat);
    Field_transient = reflect_mirror(Field_transient,obj.i_end);
    Field_transient = Propagate_E(Field_transient,obj.propagation_mat);
    
    Field_transient = Field_transient * obj.resonance_phase;
    Field_transient = reflect_mirror(Field_transient,obj.i_input);
end

obj.Field_circ = Field_total;

%------------------------------------------------------------------
% Calculate the transmitted and reflected field

Field_temp = Propagate_E(Field_total,obj.propagation_mat);
obj.Field_trans = Transmit_Reflect_Optic(Field_temp,obj.i_end);

Field_temp = Propagate_E(Field_total,obj.propagation_mat);
Field_temp = reflect_mirror(Field_temp,obj.i_end);
Field_temp = Propagate_E(Field_temp,obj.propagation_mat);
Field_temp = Field_temp * obj.resonance_phase;

Field_temp = Transmit_Reflect_Optic(Field_temp,obj.i_input);

obj.Field_ref = Field_reflec + Field_temp;

% Go back outside the substrate, usually in vacuum

if isa(obj.i_input, 'Interface')
    obj.Field_ref =  Change_E_n(obj.Field_ref,obj.i_input.n1);
end

if isa(obj.i_end, 'Interface')
    obj.Field_trans =  Change_E_n(obj.Field_trans,obj.i_end.n1);
end

%-------------------------------------------------------------------
% Calculate the round trip loss
% First get the HR interface
if isa(obj.i_input,'Mirror')
    obj.i_input = obj.i_input.I_HR;
end

if isa(obj.i_end,'Mirror')
    obj.i_end = obj.i_end.I_HR;
end

field_tmp = Field_total;
field_tmp = Normalise_E(field_tmp);
field_tmp = Propagate_E(field_tmp,obj.propagation_mat);
field_tmp = reflect_mirror(field_tmp,obj.i_end,'Ref',1);
field_tmp = Propagate_E(field_tmp,obj.propagation_mat);
field_tmp = reflect_mirror(field_tmp,obj.i_input,'Ref',1);

obj.loss_rtl =  (1 - calculate_power(field_tmp));

end