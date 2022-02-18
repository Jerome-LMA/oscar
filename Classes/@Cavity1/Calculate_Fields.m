function varargout = Calculate_Fields(Cin,varargin)
% Cout = Calculate_fields(Cin) calculate the circulating, reflected and transmitted fields
% Function used to calculated the fields inside the cavity. The laser beam must be defined outside the cavity in order to calculate the reflected field.

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
    Accuracy = 0.0001;
end
% Calculate the number of iteration to reach the steady state


RT_loss = Cin.I_input.r*Cin.I_end.r;
% Have to solve RT_loss^num_iter < 0.5*accuracy
num_iter = log(0.5*Accuracy)/(log(RT_loss));

if ~isempty(p.Results.iter) % overide the numbers of iteration
    num_iter = p.Results.iter;
end
num_iter = round(num_iter);

% The laser starts outside the input mirror, change n from 1 to mirror
% substrate refractive index

if isa(Cin.I_input, 'Interface')
    Field_in =  Change_E_n(Cin.Laser_in,Cin.I_input.n2);
else
    Field_in =  Cin.Laser_in;
end

[Field_transient,Field_reflec] = Transmit_Reflect_Optic(Field_in,Cin.I_input,'AR');
Field_total = Normalise_E(Field_transient,'Power',0);

Power_buildup = zeros(1,num_iter);

for q = 1:num_iter
    Power_buildup(q) = Calculate_Power(Field_total);
    Field_total = Field_total + Field_transient;
    Field_transient = Propagate_E(Field_transient,Cin.Propagation_mat);
    Field_transient = Reflect_Mirror(Field_transient,Cin.I_end);
    Field_transient = Propagate_E(Field_transient,Cin.Propagation_mat);
    
    Field_transient = Field_transient * Cin.Resonance_phase;
    Field_transient = Reflect_Mirror(Field_transient,Cin.I_input);
end

Cout.Field_circ = Field_total;

%------------------------------------------------------------------
% Calculate the transmitted and reflected field

Field_temp = Propagate_E(Field_total,Cin.Propagation_mat);
Cout.Field_trans = Transmit_Reflect_Optic(Field_temp,Cin.I_end);

Field_temp = Propagate_E(Field_total,Cin.Propagation_mat);
Field_temp = Reflect_Mirror(Field_temp,Cin.I_end);
Field_temp = Propagate_E(Field_temp,Cin.Propagation_mat);
Field_temp = Field_temp * Cin.Resonance_phase;

Field_temp = Transmit_Reflect_Optic(Field_temp,Cin.I_input);

Cout.Field_ref = Field_reflec + Field_temp;

% Go back outside the substrate, usually in vacuum

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

field_tmp = Field_total;
field_tmp = Normalise_E(field_tmp);
field_tmp = Propagate_E(field_tmp,Cin.Propagation_mat);
field_tmp = Reflect_Mirror(field_tmp,Cin.I_end,'Ref',1);
field_tmp = Propagate_E(field_tmp,Cin.Propagation_mat);
field_tmp = Reflect_Mirror(field_tmp,Cin.I_input,'Ref',1);

Cout.Loss_RTL =  (1 - Calculate_Power(field_tmp));


switch nargout
    case 0
    case 1
        varargout{1} = Cout;
    case 2
        varargout{1} = Cout;
        varargout{2} = Power_buildup;
    otherwise
        error('Get_info(): Too many output argument')
end

end