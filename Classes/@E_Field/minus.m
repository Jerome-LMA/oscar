function [Eout] = 	minus(E1,E2)
%minus: Overload the - function for E_Field
%simply substract the complex field of the input together

if (E1. Refractive_index ~= E2. Refractive_index)
    error('minus(): error the 2 inputs fields are taken in different media')
end

if (E1. Wavelength ~= E2. Wavelength)
    error('minus(): error the 2 inputs fields are taken in different media')
end

if (E1. Frequency_Offset ~= E2. Frequency_Offset)
    error('minus(): error the 2 inputs fields have different frequency offset')
end

if (E1.Grid ~= E2. Grid)
    error('minus(): error the 2 inputs fields are defined on different grid')
end

Eout = E1;

if isempty(E1.Field_SBl)
    
    Eout.Field = E1.Field -  E2.Field;
    
else
    
    Eout.Field = E1.Field -  E2.Field;
    Eout.Field_SBl = E1.Field_SBl - E2.Field_SBl;
    Eout.Field_SBu = E1.Field_SBu - E2.Field_SBu;
    
end

end

