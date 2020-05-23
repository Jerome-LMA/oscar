function [result, coeff] = ZernikeDecomposition(mat, degree)
% [result, coeff] = ZernikeDecomposition(mat, degree)
% Performs the decomposition of the square matrix `mat'
% into Zernike polynomials up to `degree'
%
% `result' is the approximation of `mat' as a sum of Zernike
% polynomials of degree from 0 to `degree'
% 
% `coeff' is a matrix containing row-wise the coefficients of 
% the decomposition into Zernike polynomials, namely:
%   0,0
%   1,-1    1,1
%   2,-2    2,0     2,2
%   3,-3    3,-1    3,1     3,3
%   ...
% Written by Massimo Galimberti - LMA - 2010


% error is mat is not square
[m_row,n_col] = size(mat);
if m_row ~= n_col 
    error('Matrix must be square'); 
end

% force degree to be non-negative
degree = abs(degree(1));

width = m_row;

result = zeros(size(mat));
coeff = NaN * ones(degree+1);

for n=0:degree
    col = 1;
    for m=-n:2:n
        C = ZernikeCoefficient(mat, n, m);
        result = result + C * ZernikePolynomial(n, m, width);
        coeff(n+1,col) = C;
        col = col + 1;
    end
end
