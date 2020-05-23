function E_out = Resample_E(E_in,G_out)
% Resample_E(E_in,G_out ) Resample the input E_Field on a different grid
% E_out = Resample_E(E_in,G_out) resample E_in on the grid G_out. The grids
% of E_in and G_out can have different number of points but must have the same length

p = inputParser;

% Check if the first argument is a E_Field
p.addRequired('E_in', @(x)isa(x, 'E_Field'));
% Check if the second argument is a Grid
p.addRequired('G_out', @(x)isa(x, 'Grid'));

p.parse(E_in,G_out);

% Create the output field with dummy value
E_out =  E_in;
E_out.Grid = G_out;


% Check if the grid coming with E_in and G_out have the same lengh
if E_in.Grid.Length ~= G_out.Length
    error('Resample_E(), the 2 input grids are different in length')
end

if E_in.Grid.Num_point > G_out.Num_point
    E_out.Field = interp2(E_in.Grid.D2_X,E_in.Grid.D2_Y,E_in.Field,G_out.D2_X,G_out.D2_Y,'nearest');
else
    E_out.Field = interp2(E_in.Grid.D2_X,E_in.Grid.D2_Y,E_in.Field,G_out.D2_X,G_out.D2_Y,'linear');
end

E_out.Field(isnan(E_out.Field)) = 0;

end

