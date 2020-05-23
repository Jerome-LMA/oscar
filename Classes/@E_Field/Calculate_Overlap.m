function varargout = Calculate_Overlap(varargin)
% Calculate the overlap integral between 2 E_field, it returns a complex number
% overlap = Calculate_Overlap(E1,E2), calculate the overlap for the carrier
% field E1, E2

switch nargin
    case 0
        error('Calculate_Overlap(): not enough arguments, at least one E_field must be given')
    case 1
        E1 = Normalise_E(varargin{1});
        if isempty(E1.Field_SBl)
            error('Calculate_Overlap(): no SB are present')
        end
        
        E1_temp = E1;
        E2_temp = E1;
        
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
                error('Calculate_Overlap(): Wrong number of output argument, should be none or 2')
        end
           
    case 2
        if isa(varargin{1}, 'E_Field') && isa(varargin{2}, 'E_Field')
            E1 = Normalise_E(varargin{1});
            E2 = Normalise_E(varargin{2});
            over_temp = sum(E1.Field .* conj(E2.Field), 'all') * (E1.Grid.Step)^2;
            
            switch nargout
                case 0
                    fprintf(' Power overlap between fields %g  \n ',abs(over_temp).^2)
                case 1
                    varargout{1} = over_temp;
                otherwise
                    error('Calculate_Overlap(): Too many output argument')
            end
        else
            error('Calculate_Overlap(): The two first arguments must be E_field')
        end
        
    otherwise
        error('Calculate_Overlap(): Invalid number of input arguments, no overlap calculation is made')
end


