function Eout = times(Ein,mat)
% times(E1,mat): Overload the .* function for E_Field with a matrix
% will multiply both the carrier and the sidebands
% the E_field must be the first argument

p = inputParser;

% Check if the first argument is an E_Field
p.addRequired('Ein', @(x)isa(x, 'E_Field'));

% Check if the second argument is a matrix
p.addRequired('mat', @(x) ismatrix(x));

p.parse(Ein,mat)

% The user has to check that the size of the matrices are compatible

Eout = Ein;
Eout.Field = Ein.Field .* mat;

if Ein.Nb_Pair_SB % if SB are present
    for ii=1:Ein.Nb_Pair_SB
        Eout.SB(ii).Field_lower = Ein.SB(ii).Field_lower .* mat;
        Eout.SB(ii).Field_upper = Ein.SB(ii).Field_upper .* mat;
    end
end

end


