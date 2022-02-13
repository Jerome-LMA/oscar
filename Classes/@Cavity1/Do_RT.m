function Eout = do_rt(obj,Ein,varargin)
% Eout = do_rt(Cin,Ein), do a round trip for the field Ein in the cavity
% Cin

p  = inputParser;
p.FunctionName = 'Do a round trip in the cavity';

% Check if the first argument is a Cavity1 instance
p.addRequired('Cin', @(x)isa(x, 'Cavity1'));

% Check if the second argument is a E_Field instance
p.addRequired('Ein', @(x)isa(x, 'E_Field'));

p.parse(obj,Ein,varargin{:})

Eout = Propagate_E(Ein,obj.propagation_mat);
Eout = reflect_mirror(Eout,obj.i_end);
Eout = Propagate_E(Eout,obj.propagation_mat);

Eout = Eout * obj.resonance_phase;
Eout = reflect_mirror(Eout,obj.i_input);


end