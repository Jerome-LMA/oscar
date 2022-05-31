function What_Is_Possible(Gin,varargin)
%What_is_possible: Check what you can do with the given grid

p  = inputParser;
p.FunctionName = 'Check what you can do with a certain grid size';

% Check if the first argument is an interface
p.addRequired('Gin', @(x)isa(x, 'Grid'));

% Check if the resolution of the grid if given
p.addParameter('Diam',[],@(x)isnumeric(x) && x>0);

% Check if the resolution of the grid if given
p.addParameter('Wavelength',[],@(x)isnumeric(x) && x>0);

p.parse(Gin,varargin{:})

if isempty(p.Results.wavelength)
    Wavelength = 1064E-9;
else
    Wavelength = p.Results.wavelength;
end

if isempty(p.Results.diam)
    Mirror_size = Gin.Length *0.8;
else
    Mirror_size = p.Results.diam;
end

fprintf('Wavelength: %g [nm] \n',Wavelength*1E9)
fprintf('Maximum angle for the beam allowed: %g [mrad] \n', Wavelength*Gin.Num_point/(2*Gin.Length)*1E3)
fprintf('Minimum radius of curvature for a mirror of diameter %g cm: %g [m] \n', Mirror_size,Mirror_size*1E-2*Gin.Step/Wavelength)

% Calculate the mininal distance for the use of the digital
% integration.
% Add a safety factor of 6
dis_min = Gin.Length^2 / (2 * Gin.Num_point * Wavelength);

fprintf('Minimal distance of propagation for digital integration [m]: %g  [m] \n', 6*dis_min)

end

