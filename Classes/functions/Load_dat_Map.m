function [Iout,Gout] = Load_dat_Map(file_name,varargin)
% Load_dat_map() read directly a *.dat file and return an Interface object
% with only little interpolation (the pixel size is given by the measurement)

p  = inputParser;
p.FunctionName = 'Load a map and create a grid';

% Check if the first argument is a string
p.addRequired('file_name' , @(x)ischar(x)  );

% Check if we should remove the tilt over a certain diameter (in meter)
p.addParameter('remove_tilt', [] , @(x)isnumeric(x) && x>0);

% Check if we should remove the tilt and curvarture over a certain diameter (in meter)
p.addParameter('remove_tilt_focus', [] , @(x)isnumeric(x) && x>0);

% Check if we should rotate the map by 90 deg and how many time, integer
% number 1 = 90 deg
p.addParameter('rotate', 0 , @(x)isnumeric(x));

% Check if we should offset the map, 2D vector in X and Y in m
p.addParameter('shift',[0 0], @(x)isnumeric(x));

% Check if there is a scaling factor for the map
p.addParameter('scale',[],@(x)isnumeric(x));

p.parse(file_name,varargin{:})

%p.Results

if ~exist(file_name, 'file') % Check if the file exist
    error(' Load_dat_map(): file %s does not exist! \n', file_name);
end

[Map_loaded, dx] = ReadZygoBinary(file_name);   % Load the map
Map_loaded = Make_square(Map_loaded);           % Make the map square if necessary

[m,~] = size(Map_loaded);

% Change the scale
if  ~isempty(p.Results.scale)
    Map_loaded = Map_loaded * p.Results.scale;
end

% To create a new grid, the number of points must be even
if rem(m,2)
    m = m+1;
end

Gout = Grid(m,m*dx);

% Create a dummy flat mirror
Iout = Interface(Gout,'RoC',Inf,'CA',m*dx);

% Add the loaded map to the flat mirror

if  ~isempty(p.Results.remove_tilt)
        Iout = Add_Map(Iout,Map_loaded,'reso',dx,'remove_tilt',p.Results.remove_tilt,'rotate',p.Results.rotate,...
            'shift',p.Results.shift,'scale',-1);    
elseif ~isempty(p.Results.remove_tilt_focus)
        Iout = Add_Map(Iout,Map_loaded,'reso',dx,'remove_tilt_focus',p.Results.remove_tilt_focus,'rotate',p.Results.rotate,...
            'shift',p.Results.shift,'scale',-1); 
else
        Iout = Add_Map(Iout,Map_loaded,'reso',dx,'rotate',p.Results.rotate,'scale',-1);
end

end

