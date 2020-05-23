function varargout = Calculate_power(Ein,varargin)
% Calculate the power in W of a Efield
% Use: Calculate_power(E_field) or Calculate_power(E_field,'include','all')

p  = inputParser;
p.FunctionName = 'Display an E_Field object';

% Check if the first argument is an E_Field
p.addRequired('Ein', @(x)isa(x, 'E_Field'));

% Check what are we want calculate
p.addOptional('include','carrier', @(x)strcmpi(x,'carrier') | ...
    strcmpi(x,'all') | strcmpi(x,'SB') );

p.parse(Ein,varargin{:})

if strcmp(p.Results.include,'carrier')
    power_temp_Car = sum(abs(Ein.Field).^2, 'all') * (Ein.Grid.Step)^2;
    power_temp_SB = 0;
elseif strcmp(p.Results.include,'all')
    power_temp_Car = sum(abs(Ein.Field).^2, 'all')* (Ein.Grid.Step)^2;
    power_temp_SB = ( sum(abs(Ein.Field_SBu).^2, 'all') + sum(abs(Ein.Field_SBl).^2, 'all') ) * (Ein.Grid.Step)^2;
elseif strcmp(p.Results.include,'SB')
    power_temp_SB =  ( sum(abs(Ein.Field_SBu).^2, 'all') + sum(abs(Ein.Field_SBl).^2, 'all') )  * (Ein.Grid.Step)^2;
    power_temp_Car = 0;
end

switch nargout
    case 0
        str = ['Power in the field(s): ' inputname(1) ' [W]:   %g  \n'];
        if strcmp(p.Results.include,'carrier')
            fprintf(str,power_temp_Car)
        elseif strcmp(p.Results.include,'all')
            fprintf(str,power_temp_Car+power_temp_SB)
        elseif strcmp(p.Results.include,'SB')
            fprintf(str,power_temp_SB)
        end
    case 1
        if strcmp(p.Results.include,'carrier')
            varargout{1} = power_temp_Car;
        elseif strcmp(p.Results.include,'all')
           varargout{1} = power_temp_Car+power_temp_SB;
        elseif strcmp(p.Results.include,'SB')
           varargout{1} = power_temp_SB;
        end
    case 2
        varargout{1} = power_temp_Car;
        varargout{2} = power_temp_SB;
    otherwise
        error('Calculate_power(): Too many output argument')
end


end


