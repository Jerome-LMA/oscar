function varargout = Check_Position_Tilt(Ein,varargin)
% Check_pos_tilt(): return the centroid of the beam and the tilt angle
% for both x and y direction
% Calculation done on the carrier

p  = inputParser;

% Check if the first argument is an E_Field
p.addRequired('Ein', @(x)isa(x, 'E_Field'));

p.parse(Ein,varargin{:})

Propa_length = Ein.Grid.Length*50;

Vec_cent(1,1) = sum(sum(Ein.Grid.D2_X.*abs(Ein.Field))) / sum(sum(abs(Ein.Field)));
Vec_cent(1,2)= sum(sum(Ein.Grid.D2_Y.*abs(Ein.Field))) / sum(sum(abs(Ein.Field)));

Eout = Propagate_E(Ein,Propa_length);

Vec_cent(2,1) = sum(sum(Eout.Grid.D2_X.*abs(Eout.Field))) / sum(sum(abs(Eout.Field)));
Vec_cent(2,2)= sum(sum(Eout.Grid.D2_Y.*abs(Eout.Field))) / sum(sum(abs(Eout.Field)));

beam_pos(1) =  Vec_cent(1,1);
beam_pos(2) =  Vec_cent(1,2);

beam_tilt(1) = (Vec_cent(2,1)-Vec_cent(1,1))/Propa_length;
beam_tilt(2) = (Vec_cent(2,2)-Vec_cent(1,2))/Propa_length;

switch nargout
    case 0
        fprintf('Beam position x [m]: %g  \t \t Beam position y [m] %g  \n',beam_pos(1),beam_pos(2))
        fprintf('Beam angle x [mrad]: %g  \t \t Beam angle y [mrad] %g  \n',beam_tilt(1)*1E3,beam_tilt(2)*1E3)
    case 1
        varargout{1} = beam_pos;
    case 2
        varargout{1} = beam_pos;
        varargout{2} = beam_tilt;
    otherwise
        error('Check_pos_tilt(): Too many output argument')
end


end

