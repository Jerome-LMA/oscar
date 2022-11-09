function Iout = Resize_Interface(Iin,G2)
% Iout = Resize_interface(In,G2) Create a new interface which Iin resampled
% on the grid G2

p  = inputParser;
p.FunctionName = 'Resampled a map on a new grid';

% Check if the first argument is an interface
p.addRequired('Iin', @(x)isa(x, 'Interface'));

% Check if the  second argument is a grid
p.addRequired('G2', @(x)isa(x, 'Grid'));

p.parse(Iin,G2)

Iout = Interface(G2);

Iout.Grid = interp2(Iin.Grid.D2_X,Iin.Grid.D2_Y,Iin.surface,G2.D2_X,G2.D2_Y);
Iout.mask = interp2(Iin.Grid.D2_X,Iin.Grid.D2_Y,Iin.mask,G2.D2_X,G2.D2_Y,'nearest');

Iout.T = Iin.T;
Iout.L = Iin.L;

Iout.n1 = Iin.n1;
Iout.n2 = Iin.n2;

Iout.r = Iin.r;
Iout.t = Iin.t;




end

