function Eout = add_tilt(Ein,tilt_angle,varargin)
% Add_tilt: Add a tilted wavefront to an E_field
% E2 = Add_tilt(E1,tilt_angle)   tilt the wavefront by tilt_angle in the
% horizontal direction, tilt is given in radian and is a small number
% E2 = Add_tilt(E1,tilt_angle,'dir','y') tilt the wavefront in the vertical
% direction

p = inputParser;
p.FunctionName = 'Add tilt on the wavefront';

% Check if the first argument is an E_Field
p.addRequired('Ein', @(x)isa(x, 'E_Field'));

% Check if the second argument is an angle
p.addRequired('tilt_angle', @(x) isnumeric(x));

% Optional parameter 'x' or 'y' to specify the axis of the tilt
p.addParameter('dir','x', @(x)strcmpi(x,'x') | strcmpi(x,'y'));

p.parse(Ein,tilt_angle,varargin{:})

% Calculate the additional tilted wavefront
if strcmp(p.Results.dir,'x')
    Add_WF = exp(-1i*Ein.k_prop.*Ein.Grid.D2_X*tilt_angle);
else
    Add_WF = exp(-1i*Ein.k_prop.*Ein.Grid.D2_Y*tilt_angle);
end

Eout = Ein.*Add_WF;
        
end

