function Cout = Optimal_Matching(Cin,varargin)
% Optimal_matching()  set the optimal mode matching for a cavity, so change
% the input beam

p  = inputParser;

% Check if the first argument is an interface
p.addRequired('Cin', @(x)isa(x, 'Cavity1'));

p.parse(Cin,varargin{:})

Cout = Cin;
beam_size_temp = Check_stability(Cin,'Display', false);

New_E_input = E_Field(Cin.Laser_in.Grid,'w',beam_size_temp(1),'R',beam_size_temp(2));

Cout.Laser_in = New_E_input;


end

