function Eout = Do_RT(Cin,Ein,varargin)
% Eout = Do_RT(Cin,Ein), do a round trip for the field Ein in the cavity
% Cin

p  = inputParser;
p.FunctionName = 'Do a round trip in the cavity';

% Check if the first argument is a Cavity1 instance
p.addRequired('Cin', @(x)isa(x, 'Cavity1'));

% Check if the second argument is a E_Field instance
p.addRequired('Ein', @(x)isa(x, 'E_Field'));

p.parse(Cin,Ein,varargin{:})

Eout = Propagate_E(Ein,Cin.Propagation_mat);
Eout = Reflect_mirror(Eout,Cin.I_end);
Eout = Propagate_E(Eout,Cin.Propagation_mat);

Eout = Eout * Cin.Resonance_phase;
Eout = Reflect_mirror(Eout,Cin.I_input);


end