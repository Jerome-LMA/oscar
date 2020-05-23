function coeff = ZernikeCoefficient(mat, n, m)
% coeff = ZernikeCoefficient(mat, n, m, normalize=false)
% Returns the coefficient of order n,m of the decomposition
% of the square matrix `mat' into Zernike polynomials
% Written by Massimo Galimberti - LMA - 2010


% error if mat is not square
[m_row,n_col] = size(mat); if m_row ~= n_col error('Matrix must be square'); end

% force n to be a non-negative integer scalar
n = floor(abs(n(1)));
% force m to be a scalar integer
m = fix(m(1));
% force m to be -n <= m <= n
if (abs(m) > n )
    m = n * sign(m);
end

% change NaN to zeros
mat( isnan(mat) ) = 0;

% Zernike polynomial
width = m_row;
Znm = ZernikePolynomial(n, m, width);

% to get the integration coefficient, we impose that the projection
% of Znm over itself be 1
int_coeff = sum( Znm .* Znm, 'all' );

% projection of the input matrix on the polynomial
coeff = sum( mat .* Znm, 'all' ) / int_coeff;
