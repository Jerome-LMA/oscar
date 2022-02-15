function Iout = Tilt_interface(Iin,AoI,varargin)
%     Tilt_interface() Project a surface according to the angle of
%     incidence of the light beam
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Iout =  Tilt_interface(Iin,AoI) Tilt the surface Iin according to the
%    angle of incidence AoI.
%    !! Only for near flat surface. In case, a curvature has to be included,
%    use the definition of the class Interface, that will include the
%    proper astigmatic curvature as well as the elliptical aperture.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p  = inputParser;
p.FunctionName = 'Tilt a map';

% Check if the first argument is an interface
p.addRequired('Iin', @(x)isa(x, 'Interface'));

% Check if the second argument is a file to load or a matrix
p.addRequired('AoI', @(x)isnumeric(x) && x>=0);

p.parse(Iin,AoI,varargin{:})

%p.Results

Iout = Iin;

Angle_tilt = AoI * pi / 180; % Pass in radian

G1 = Iin.Grid;

% Create a streched grid
[NG_D2_X,NG_D2_Y] = meshgrid(G1.Axis*cos(Angle_tilt),G1.Axis);

Iout.surface = interp2(G1.D2_X,G1.D2_Y,Iin.surface,NG_D2_X,NG_D2_Y,'spline',0);

end
