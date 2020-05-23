function varargout = Calculate_Overlap_SB_Car(Ein,varargin)
% Calculate the overlap integral between the SB and the carrier, it returns
% two complex number one per SB
% overlap = Calculate_Overlap(E1), calculate the overlap between the carrier
% and the SB fields

p = inputParser;
p.FunctionName = 'Display an E_Field object';

% Check if the first argument is an E_Field
p.addRequired('Ein', @(x)isa(x, 'E_Field'));

p.parse(Ein,varargin{:})

E1_temp = Ein;
E2_temp = Ein;

E1_temp.Field = E1_temp.Field;
E2_temp.Field = E1_temp.Field_SBl;
O_lower = Calculate_Overlap(E1_temp,E2_temp);

E2_temp.Field = E1_temp.Field_SBu;
O_upper = Calculate_Overlap(E1_temp,E2_temp);

switch nargout
    case 0
        fprintf('  Power overlap between carrier and lower sideband %g  \n',abs(O_lower).^2)
        fprintf('  Power overlap between carrier and upper sideband %g  \n',abs(O_upper).^2)
    case 2
        varargout{1} = O_lower;
        varargout{2} = O_upper;
    otherwise
        error('Calculate_Overlap_SB_Car(): Wrong number of output argument, should be none or 2')
end

end


