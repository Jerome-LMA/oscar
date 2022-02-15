function Iout = Resize_interface(In,G2)
% Iout = Resize_interface(In,G2) Create a new interface which Iin resampled
% on the grid G2

p  = inputParser;
p.FunctionName = 'Resampled a map on a new grid';

% Check if the first argument is an interface
p.addRequired('Iin', @(x)isa(x, 'Interface'));

% Check if the  second argument is a grid
p.addRequired('G2', @(x)isa(x, 'Grid'));

p.parse(In,G2)

Iout = Interface(G2);





       Grid: [1x1 Grid]
    surface: [742x742 double]
       mask: [742x742 double]
          T: 0.1000
          L: 0
         n1: 1
         n2: 1.4500
          t: 0 + 0.3162i
          r: 0.9487



end

