function prop_out = Refine_DI_prop(prop_in,E_in)
%Refine_DI_prop: try to avoid the aliasing effect for the digital integration
% prop_out = Refine_DI_prop(prop_in) calculate the propagation operator for
% the digital integration on a grid size with a larger number of points.
% !! That may also reduce the accuracy of the calculations !!

p  = inputParser;
p.FunctionName = 'Modify the operator for the digital integration';

% Check if the first argument is a E_Field object
p.addRequired('prop_in', @(x)isa(x, 'Prop_operator'));
p.addRequired('E_in', @(x)isa(x, 'E_Field'));
p.parse(prop_in,E_in);

prop_out = prop_in;

% Define a new grid with more points but the same physical size
G2 = Grid(2048,prop_in.Grid.Length);

% Define the propagation distance
%Prop.dist = prop_in.dist;

%Ratio_grid = G2.Num_point / G1.Num_point;

Prop.mat_DI = zeros(2*G2.Num_point-1);

X_tmp_1D = linspace(2*G2.Axis(1),-2*G2.Axis(1),2*G2.Num_point-1);
Y_tmp_1D = linspace(2*G2.Axis(1),-2*G2.Axis(1),2*G2.Num_point-1);

[X_tmp_2D Y_tmp_2D] = meshgrid(X_tmp_1D,Y_tmp_1D);
tmp_rad = sqrt(X_tmp_2D.^2 + Y_tmp_2D.^2 + prop_in.dist^2);

Prop.mat_DI = 1/(2*pi) * (exp(-1i * E_in.k_prop *tmp_rad) ./ tmp_rad) .* (prop_in.dist ./  tmp_rad)  .* (1 ./  tmp_rad + 1i * E_in.k_prop );
Prop.mat_DI = Prop.mat_DI * (G2.Step)^2;

Prop.mat_DI = fftshift(fft2(Prop.mat_DI));

% Cut a square in the big matrix
Central_point = G2.Num_point;
X_cut = (Central_point - (prop_in.Grid.Num_point-1)):(Central_point + (prop_in.Grid.Num_point-1));
Y_cut = (Central_point - (prop_in.Grid.Num_point-1)):(Central_point + (prop_in.Grid.Num_point-1));

prop_out.mat_DI = ifftshift(Prop.mat_DI(X_cut,Y_cut));


% Comparison of the 2 matrices

% figure(11)
% imagesc(abs(prop_in.mat_DI)); axis square
% 
% figure(12)
% imagesc(abs(prop_out.mat_DI)); axis square

end

