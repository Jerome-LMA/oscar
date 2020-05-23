function E = Add_tilt(varargin)
% Add_tilt: Add a tilted wavefront to an E_field
% E2 = Add_tilt(E1,tilt_angle)   tilt the wavefront by tilt_angle in the
% horizontal direction
% E2 = Add_tilt(E1,tilt_angle,'y') tilt the wavefront in the vertical
% direction

switch nargin
    case {0,1}
        disp('Add_tilt(): not enough arguments, at least an object E_field and an angle must be given')
        return
    case 2
        E = varargin{1};
        tilt_angle =  varargin{2};
        E.Field = E.Field .* exp(-1i*E.k_prop.*E.Grid.D2_X*tilt_angle);
        
        if ~isempty(E.Field_SBl)
            E.Field_SBl = E.Field_SBl .* exp(-1i*E.k_prop.*E.Grid.D2_X*tilt_angle);
            E.Field_SBu = E.Field_SBu .* exp(-1i*E.k_prop.*E.Grid.D2_X*tilt_angle);
        end
        
    case 3
        if strcmp(varargin{3},'x')
            E = varargin{1};
            tilt_angle =  varargin{2};
            E.Field = E.Field .* exp(-1i*E.k_prop.*E.Grid.D2_X*tilt_angle);
            
            if ~isempty(E.Field_SBl)
                E.Field_SBl = E.Field_SBl .* exp(-1i*E.k_prop.*E.Grid.D2_X*tilt_angle);
                E.Field_SBu = E.Field_SBu .* exp(-1i*E.k_prop.*E.Grid.D2_X*tilt_angle);
            end
            
            
        elseif strcmp(varargin{3},'y')
            E = varargin{1};
            tilt_angle =  varargin{2};
            E.Field = E.Field .* exp(-1i*E.k_prop.*E.Grid.D2_Y*tilt_angle);
            
            if ~isempty(E.Field_SBl)
                E.Field_SBl = E.Field_SBl .* exp(-1i*E.k_prop.*E.Grid.D2_Y*tilt_angle);
                E.Field_SBu = E.Field_SBu .* exp(-1i*E.k_prop.*E.Grid.D2_Y*tilt_angle);
            end
                     
        else
            disp('Wrong third argument for the tilt it must be the string x or y')
            return
        end
        
    otherwise
        disp('Add_tilt(): Invalid number of input arguments, no tilt has been added')
        return
        
end

