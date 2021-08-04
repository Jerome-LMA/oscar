% Find the two RoCs of a linear Fabry Perot cavity for a given spot size on
% the mirrors. All in SI units.

% example for Advanced Virgo arm cavity

Wavelength = 2000E-9;    % lambda = 1064 nm
Length = 10E3;            % Cavity length = 3 km

w_IM = 90E-3;          % beam radius on the input mirror  = 48.7 mm
w_EM = 90E-3;          % beam radius on the end mirror  = 58.0 mm

tmp = Calculate_beam_size(Wavelength,Length,1421,1682.1);

fun2zero = @(x) Calculate_beam_size(Wavelength, Length, x(1), x(2)) - [w_IM,w_EM];

options = optimoptions('fsolve','FiniteDifferenceType','central','OptimalityTolerance',1E-12,'Display','iter'); %'FiniteDifferenceStepSize',0.01
[x,fval] = fsolve(fun2zero, [Length*0.8, Length*0.8],options);

fprintf('RoC IM [m]: %g \n',x(1))
fprintf('RoC EM [m]: %g \n',x(2))

g1 = 1 -  Length/x(1);
g2 = 1 -  Length/x(2);

w0 = sqrt((Wavelength * Length / pi) * sqrt ((g1*g2*(1-g1*g2))/(g1+g2-2*g1*g2)^2));
distITM_waist = (g2*(1 - g1)*Length )/(g1+g2 - 2*g1*g2);


function beam_size = Calculate_beam_size(lambda, Cavity_length, RoC_IM, RoC_EM)

g1 = 1 -  Cavity_length/RoC_IM;
g2 = 1 -  Cavity_length/RoC_EM;

%w0 = sqrt((lambda * Cavity_length / pi) * sqrt ((g1*g2*(1-g1*g2))/(g1+g2-2*g1*g2)^2));
%distITM_waist = (g2*(1 - g1)*Cavity_length )/(g1+g2 - 2*g1*g2);

w_IM = sqrt( (lambda * Cavity_length / pi) * sqrt (g2 /(g1*(1-g1*g2))));
w_EM = sqrt( (lambda * Cavity_length / pi) * sqrt (g1 /(g2*(1-g1*g2))));

beam_size(1) = w_IM;
beam_size(2) = w_EM;

end

