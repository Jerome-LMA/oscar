function coeff = LGcoefficient(mat, p, l, lambda, w, curv, gridsize)
% coeff = LGcoefficient(mat, p, l, lambda, w, curv, gridsize)
%
% Returns the complex coefficient of order p,l of the decomposition
% of the square matrix `mat' into Laguerre-Gauss modes
%
% p, l              order of the mode
% lambda            wavelength
% w                 spot size of the beam
% curv              beam curvature (1/ROC)
% gridsize          grid size


% error if mat is not square
if issquare(mat) == 0
    error('input matrix must be square');
end

% force p to be non-negative
p = abs(p);

% normalized LG mode
num_gridpoints = rows(mat);
LGmn = LaguerreGauss(p, l, lambda, w, curv, gridsize, num_gridpoints);
% size of the single pixel
dx = gridsize / num_gridpoints;

% projection of the input matrix on the polynomial
coeff = sum( mat .* conj(LGmn) .* dx^2, 'all' );
