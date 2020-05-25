function [Eout] = plus(E1,E2)
%Plus: Overload the + function for E_Field
%simply add the complex field of the input together

if (E1. Refractive_index ~= E2. Refractive_index)
    error('plus(): error the 2 inputs fields are taken in different media')
end

if (E1. Wavelength ~= E2. Wavelength)
    error('plus(): error the 2 inputs fields are taken in different media')
end

if (E1.Nb_Pair_SB ~= E2.Nb_Pair_SB)
    error('plus(): error the 2 inputs fields have different number of sidebands')
end

if (E1.Grid ~= E2. Grid)
    error('plus(): error the 2 inputs fields are defined on different grid')
end


Eout = E1;

Eout.Field = E1.Field +  E2.Field;

if Eout.Nb_Pair_SB
    for ii=1:Eout.Nb_Pair_SB  % if SB are present
        if (E1.SB(ii).Frequency_Offset ~= E2.SB(ii).Frequency_Offset)
            error('plus(): error the 2 input fields have sidebands frequency')
        end
        Eout.SB(ii).Field_lower = E1.SB(ii).Field_lower + E2.SB(ii).Field_lower;
        Eout.SB(ii).Field_upper = E1.SB(ii).Field_upper + E2.SB(ii).Field_upper;
    end
end

end

