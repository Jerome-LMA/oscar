function varargout = Calculate_Power(Ein,varargin)
% Calculate the power in W of a E_field, it could be carrier and/or
% sidebands
% Use: Calculate_power(E_field) or Calculate_power(E_field,'include','all')

p  = inputParser;

% Check if the first argument is an E_Field
p.addRequired('Ein', @(x)isa(x, 'E_Field'));

% Check what are we want calculate
p.addOptional('include','carrier', @(x)strcmpi(x,'carrier') | ...
    strcmpi(x,'all') | strcmpi(x,'SB') );

% Check the number of the sidebands we want to deal with
p.addParameter('SB_num',1, @(x) isnumeric(x)  && (x>0) && (mod(x,1) == 0)); % check if the number of the SB pair is positive and integer

p.parse(Ein,varargin{:})
SB_number = p.Results.SB_num;

% Initiate the variables:
power_temp_Car = 0;
power_temp_SB_lower = 0;
power_temp_SB_upper = 0;

% check if Ein is on the GPU, if yes do the calculations in the GPU

if exist('isgpuarray','file') % to be compatible with version < 2020b
    Run_on_GPU = isgpuarray(Ein.Field);
else
    Run_on_GPU = false;
end

if strcmp(p.Results.include,'carrier')
    if Run_on_GPU
        power_temp_Car = sum(abs(arrayfun(@times,Ein.Field,Ein.Field) ),'all');
    else
        power_temp_Car = sum(abs(Ein.Field(:)).^2);
    end
    
    power_temp_Car = power_temp_Car * Ein.Grid.Step_sq;
    power_temp_SB = 0;
    
elseif strcmp(p.Results.include,'all')
    power_temp_Car = sum(abs(Ein.Field(:).^2))* (Ein.Grid.Step)^2;
    if Ein.Nb_Pair_SB
        for ii = 1:Ein.Nb_Pair_SB
            power_temp_SB_lower = power_temp_SB_lower + sum(abs(Ein.SB(SB_number).Field_lower(:)).^2)* (Ein.Grid.Step)^2;
            power_temp_SB_upper = power_temp_SB_upper + sum(abs(Ein.SB(SB_number).Field_upper(:)).^2)* (Ein.Grid.Step)^2;
        end
        
    end
    
elseif strcmp(p.Results.include,'SB')
    if Ein.Nb_Pair_SB
        if SB_number > Ein.Nb_Pair_SB
            error('Calculate_Power(): SB number invalid')
        end
        power_temp_SB_lower =  sum(abs(Ein.SB(SB_number).Field_lower(:)).^2)* (Ein.Grid.Step)^2;
        power_temp_SB_upper =  sum(abs(Ein.SB(SB_number).Field_upper(:)).^2)* (Ein.Grid.Step)^2;
    end
    power_temp_Car = 0;
end

switch nargout
    case 0
        str = ['Power in the field(s): ' inputname(1) ' [W]:   %g  \n'];
        if strcmp(p.Results.include,'carrier')
            fprintf(str,power_temp_Car)
        elseif strcmp(p.Results.include,'all')
            fprintf(str,power_temp_Car+power_temp_SB_lower+power_temp_SB_upper)
        elseif strcmp(p.Results.include,'SB')
            fprintf(str,power_temp_SB_lower+power_temp_SB_upper)
        end
    case 1
        if strcmp(p.Results.include,'carrier')
            varargout{1} = power_temp_Car;
        elseif strcmp(p.Results.include,'all')
            varargout{1} = power_temp_Car+power_temp_SB_lower+power_temp_SB_upper;
        elseif strcmp(p.Results.include,'SB')
            varargout{1} = power_temp_SB_lower+power_temp_SB_upper;
        end
    case 2
        varargout{1} = power_temp_SB_lower;
        varargout{2} = power_temp_SB_upper;
    otherwise
        error('Calculate_power(): Too many output argument')
end


end


