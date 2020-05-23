function H = HermitePolynomial(n, x)
% H = HermitePolynomial(n, x)
%
% compute the Hermite polynomial of order n
% x is the polynomial variable (might be a scalar, a vector, or a matrix)
%Author: Massimo Galimberti 2010


% force arguments to be meaningful
n = floor(abs(n(1)));	% n must be a non-negative integer

switch n 
	case {0}
		H = 1;

	case {1}
		H = 2*x;

	case {2}
		H = 4*x.^2 - 2;

	otherwise
		H = 2*x .* HermitePolynomial(n-1,x) - 2*(n-1) * HermitePolynomial(n-2,x);
end