function Z = ZernikePolynomial(n, m, width, normalize)
% Z = ZernikePolynomial(n, m, width, normalize=false)
% Generates a square matrix containing the Zernike polynomial
% of degree n and order m on a unitary circle whose diameter
% is the matrix size
%
% n         degree of the function (n = 0,1,2,...)
% m         order of the function (m = -n,...,n)
% width     size of the output matrix
% normalize if `true', the polynomial is normalized so that
%           the integral of Z^2 equals 1
% Written by Massimo Galimberti - LMA - 2010


% force n to be non-negative
n = fix( abs( n(1) ) );
% force m to be -n <= m <= n
m = fix( m(1) );
if (abs(m) > n )
    m = n * sign(m);
end
% force width to be positive
if width == 0
    width = 1;
else
    width = fix( abs( width(1) ) );
end
% check for the `normalize' argument, otherwise default it to false
if nargin() < 4
    normalize = false;
end


% define the normalized radial and azimuthal coordinate
x = linspace(-1, 1, width);
y = x';
[xx, yy] = meshgrid(x,y);
cplx = xx + 1i*yy;
rho   = abs(cplx);
theta = angle(cplx);

% Zernike Polynomial
R = ZernikeRadialPolynomial(n, abs(m), rho);

% crop outside the unitary circle
idx = find(rho>1);
R(idx) = 0;

% output
if (m == 0)
    Z = R;
elseif (m>0)
    Z = R .* cos(m * theta);
else
    Z = R .* sin(abs(m) * theta);
end

% the polynomial is adjusted by adding a constant factor so that:
% - the integral of Z00 equals pi
% - the integral of any other Znm equals zero
if ~(n==0 && m==0)
    const = sum(Z, 'all');
    npixel = numel(Z)-numel(idx);
    Z = Z - const/npixel;
    Z(idx) = 0;
end

% normalization
if normalize
    norm_coeff = pi/2/(n+1);
    if m==0; norm_coeff=norm_coeff*2; end
    Z = Z / sqrt(norm_coeff);
end
