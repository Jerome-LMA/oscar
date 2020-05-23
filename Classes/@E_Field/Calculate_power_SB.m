function varargout = Calculate_power_SB(varargin)
% Calculate the power in W of the sidebands of an Efield
% Use Calculate_power(E_field)
switch nargin
    case 0
        error('Calculate_power_SB(): not enough arguments, at least an object E_field must be given')
    case 1
        if isa(varargin{1}, 'E_Field')
            if isempty(varargin{1}.Field_SBl)
                error('Calculate_power_SB(): no sidebands are present')
            end
            power_temp1 = sum(abs(varargin{1}.Field_SBl).^2, 'all') * (varargin{1}.Grid.Step)^2;
            power_temp2 = sum(abs(varargin{1}.Field_SBu).^2, 'all') * (varargin{1}.Grid.Step)^2;
            
            switch nargout
                case 0
                    str = ['Power in the lower sideband of the field ' inputname(1) ' [W]:   %g  \n'];
                    fprintf(str,power_temp1)
                    str = ['Power in the upper sideband of the field ' inputname(1) ' [W]:   %g  \n'];
                    fprintf(str,power_temp2)
                case 1
                    varargout{1} = power_temp1;
                case 2
                    varargout{1} = power_temp1;
                    varargout{2} = power_temp2;
                otherwise
                    error('Calculate_power(): Too many output argument')
            end
            
        else
            error('Calculate_power(): The first argument must be an Efield')
        end
    otherwise
        error('Calculate_power(): Invalid number of input arguments, no power calculation is made')
end


