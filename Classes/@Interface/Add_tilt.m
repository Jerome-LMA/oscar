function I = Add_tilt(varargin)
% Add_tilt: Add a tilted wavefront to an E_field (or interface)
% The tilt angle in radian
% E2 = Add_tilt(E1,tilt_angle) tilt the wavefront (or surface) by tilt_angle in the
% horizontal direction
% E2 = Add_tilt(E1,tilt_angle,'y') tilt the wavefront (or surface) in the vertical
% direction

switch nargin
    case {0,1}
        disp('Add_tilt(): not enough arguments, at least an object E_field (or surface) and an angle must be given')
        return
    case 2
        I = varargin{1};
        tilt_angle =  varargin{2};
        I.surface = I.surface + I.Grid.D2_X*sin(tilt_angle);
    case 3
        if strcmp(varargin{3},'x')
            I = varargin{1};
            tilt_angle =  varargin{2};
            I.surface = I.surface + I.Grid.D2_X*sin(tilt_angle);
            
        elseif strcmp(varargin{3},'y')
            I = varargin{1};
            tilt_angle =  varargin{2};
            I.surface = I.surface + I.Grid.D2_Y*sin(tilt_angle);
        else
            disp('Wrong third argument for the tilt it must be the string x or y')
            return
        end
              
    otherwise
        disp('Add_tilt(): Invalid number of input arguments, no tilt has been added')
        return
        
end

