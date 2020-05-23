function L = LaguerrePolynomial(p, l, x)
% compute the generalized Laguerre polynomial of order p,l
% x is the polynomial variable (might be a scalar, a vector, or a matrix)
% Author: Massimo Galimberti 2010


% force arguments to be meaningful
p = floor(abs(p(1)));	% p must be a non-negative integer
l = floor(l(1));		% l must be an integer

switch p 
	case {0}
		L = 1;

	case {1}
		L = l + 1 - x;

	case {2}
		L = 0.5*(l+1)*(l+2) - (l+2)*x + 0.5*x.^2;

	otherwise
		L = 1/p * ( (2*p+l-1-x).*LaguerrePolynomial(p-1,l,x) - (p+l-1).*LaguerrePolynomial(p-2,l,x) );
end