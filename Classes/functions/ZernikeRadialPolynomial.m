function R = ZernikeRadialPolynomial(n, m, rho)
% R = ZernikePolynomial(n, m, rho, norm=false)
% Generates a matrix the same size of rho containing
% the Zernike radial polynomial of degree n and order m
% Returns a zero matrix if n-m is odd
%
% n         degree of the function (n = 0,1,2,...)
% m         order of the function (m = 0,...,n)
% rho       radial coordinate (may be a vector or a matrix)
% Written by Massimo Galimberti - LMA - 2010


% force n to be non-negative
n = fix( abs( n(1) ) );
% force m to be non-negative and <= n
m = min( fix (abs ( m(1) ) ) , n );

R = zeros( size(rho) );

% return zero if n-m is odd
if mod(m-n, 2) == 0
    for k = 0 : ((n-m)/2);
        coeff = (-1).^k .* factorial(n-k) ./ factorial(k) ./ factorial((n+m)/2-k) ./ factorial((n-m)/2-k);
        R = R + coeff .* rho.^(n-2*k);
    end
end
