function [Eout] = mtimes(E1,E2)
%mtimes(E1,E2): Overload the * function for E_Fields
% use: to multiply 2 E_field or one E_field and a scalar (or a scalar with
% an E_field)

if isa(E1, 'E_Field') % First E1 is a E_field
    
    if isa(E2, 'E_Field')
        
        if (E1. Refractive_index ~= E2. Refractive_index)
            disp('mtimes(): error the 2 inputs fields are taken in different media')
            return
        end
        
        if (E1. Wavelength ~= E2. Wavelength)
            disp('mtimes(): error the 2 inputs fields are taken in different media')
            return
        end
        
        if (E1. Frequency_Offset ~= E2. Frequency_Offset)
            disp('mtimes(): error the 2 inputs fields have different frequency offset')
            return
        end
        
        if (E1.Grid ~= E2. Grid)
            disp('mtimes(): error the 2 inputs fields are defined on different grid')
            return
        end
        
        Eout = E1;
        Eout.Field = E1.Field .* E2.Field;
        
    elseif isscalar(E2)
        
        if isempty(E1.Field_SBl)
            Eout = E1;
            Eout.Field = E1.Field * E2;
            
        else
            Eout = E1;
            Eout.Field = E1.Field * E2;
            Eout.Field_SBl = E1.Field_SBl * E2;
            Eout.Field_SBu = E1.Field_SBu * E2;
            
        end
    else
        error('mtimes(): The second input must be either a E_Field or a scalar')
    end
    
elseif isscalar(E1)
    
    if isa(E2, 'E_Field')
        
        if isempty(E2.Field_SBl)
            Eout = E2;
            Eout.Field = E2.Field * E1;
            
        else
            Eout = E2;
            Eout.Field = E2.Field * E1;
            Eout.Field_SBl = E2.Field_SBl * E1;
            Eout.Field_SBu = E2.Field_SBu * E1;
            
        end
    else
        error('mtimes(): The second input must be either a E_Field or a scalar')
    end
    
else
    error('mtimes(): The first input must be either a E_Field or a scalar')
end

end


