function varargout = Calculate_Overlap(Ein,varargin)
% overlap = Calculate_Overlap(E1,E2), calculate the overlap for the carrier
% between the fields E1, E2
% Calculate the overlap integral between 2 E_field, it returns a complex number
% without output, calculate the power overlap in the command window
% if one input, calculate the overlap between carrier and sidebands

p = inputParser;

% Check if the first argument is an E_Field
p.addRequired('Ein', @(x)isa(x, 'E_Field'));

% Check if the second argument is an E_Field
p.addOptional('E2',[],@(x)isa(x, 'E_Field') || isempty(x));

% Check the number of the sidebands we want to deal with
p.addParameter('SB_num',1, @(x) isnumeric(x)  && (x>0) && (mod(x,1) == 0)); % check if the number of the SB pair is positive and integer

p.parse(Ein,varargin{:})
SB_number = p.Results.SB_num;
E2 = p.Results.E2;

% Proper normalisation for the beam power from the electric field
% Power in W = 1/2 * c * epsilon 0 * E^2
Cte_conversion = (1/2) * 3.00E8 * 8.85E-12;


if exist('isgpuarray','file') % to be compatible with version < 2020b
    Run_on_GPU = isgpuarray(Ein.Field);
else
    Run_on_GPU = false;
end

if isempty(E2) % in that case do the overlap between the carrier and SBs
    
    if ~Ein.Nb_Pair_SB
        error('Calculate_Overlap(): no sidebands are present, it needs a second argument')
    end
    
    E2_temp = Ein;
    
    E2_temp.Field = Ein.SB(SB_number).Field_lower;
    O_lower = Calculate_Overlap(Ein,E2_temp);
    
    E2_temp.Field = Ein.SB(SB_number).Field_upper;
    O_upper = Calculate_Overlap(Ein,E2_temp);
    
    switch nargout
        case 0
            fprintf('  Power overlap between carrier and lower sideband: %g  \n',abs(O_lower).^2)
            fprintf('  Power overlap between carrier and upper sideband: %g  \n',abs(O_upper).^2)
        case 1
            varargout{1} = 0.5*(O_lower +  O_upper);
        case 2
            varargout{1} = O_lower;
            varargout{2} = O_upper;
        otherwise
            error('Calculate_Overlap(): Wrong number of output argument, should be none or 2')
    end
    
    
else
    Ein = Normalise_E(Ein);
    E2 = Normalise_E(E2);
    
    if Run_on_GPU
        over_temp = sum(arrayfun(@times,Ein.Field,conj(E2.Field)),'all') * Ein.Grid.Step_sq;
    else
        over_temp = sum(Ein.Field .* conj(E2.Field),'all') * Ein.Grid.Step_sq;
    end
    
    over_temp = over_temp*Cte_conversion;
    
    switch nargout
        case 0
            fprintf('  Power overlap between fields: %g  \n ',abs(over_temp).^2)
        case 1
            varargout{1} = over_temp;
        otherwise
            error('Calculate_Overlap(): Too many output argument')
    end
end

end


