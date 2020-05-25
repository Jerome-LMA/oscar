function Eout = Change_E_n(Ein,new_n)
% Change_E_nEin,new_n): change the refractive index of the
% field. Simple formula
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change_E_refractive_index(Ein,new_n) will change the refractive index of
% the E_field, equivalent to a flat interface between two media
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Ein.Refractive_index == new_n
    error('Change_E_n(): old and new refractive index are the same ')
end

Eout = Ein;

% Update the output
Eout.Refractive_index = new_n;
Eout.k_prop = Ein.k_prop * (new_n/Ein.Refractive_index);

end