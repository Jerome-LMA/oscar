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
                
        if (E1.Grid ~= E2. Grid)
            disp('mtimes(): error the 2 inputs fields are defined on different grid')
            return
        end
        
        %Eout = E1;
        %Eout.Field = E1.Field .* E2.Field;
        
        error('mtimes(): multiplication between 2 fields not supported since too ambiguous')
        
    elseif isscalar(E2)
        Eout = E1;
        Eout.Field = E1.Field * E2;
        
        if E1.Nb_Pair_SB % if SB are present
            for ii=1:E1.Nb_Pair_SB
                Eout.SB(ii).Field_lower = E1.SB(ii).Field_lower * E2;
                Eout.SB(ii).Field_upper = E1.SB(ii).Field_upper * E2;
            end
        end
        
    else
        error('mtimes(): The second input must be either a E_Field or a scalar')
    end
    
elseif isscalar(E1)
    
    if isa(E2, 'E_Field')
        Eout = E2;
        Eout.Field = E2.Field * E1;
        if E2.Nb_Pair_SB % if SB are present
            for ii=1:E2.Nb_Pair_SB
                Eout.SB(ii).Field_lower = E2.SB(ii).Field_lower * E1;
                Eout.SB(ii).Field_upper = E2.SB(ii).Field_upper * E1;
            end
        end
        
    else
        error('mtimes(): The second input must be either a E_Field or a scalar')
    end
    
else
    error('mtimes(): The first input must be either a E_Field or a scalar')
end

end


