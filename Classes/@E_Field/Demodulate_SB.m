function [p,q] = Demodulate_SB(Ein,varargin)
% [p q] = Demodulate_SB(E1), demodulate the carrier with the sidebdands
% [p q] are the signal in phase and in quadrature

p = inputParser;
p.FunctionName = 'Demodulate a signal';

% Check if the first argument is a E_Field
p.addRequired('Ein', @(x)isa(x, 'E_Field'));

% Check the number of the sidebands we want to deal with
p.addParameter('SB_num',1, @(x) isnumeric(x)  && (x>0) && (mod(x,1) == 0)); % check if the number of the SB pair is positive and integer

% Check if a demodulation phase is given (must be in radian)
p.addParameter('phase',0,@(x)isnumeric(x));

p.parse(Ein,varargin{:})

if ~Ein.Nb_Pair_SB
    error('Demodulate_SB(): no sidebands are present')
end

SB_number = p.Results.SB_num;

% Check the number of SB fields to be correct
if SB_number > Ein.Nb_Pair_SB
    error('Demodulate_SB(): requested SB field not present')
end

Esignal = sum(sum( Ein.Field .* conj(Ein.SB(SB_number).Field_upper) + conj(Ein.Field) .* Ein.SB(SB_number).Field_lower)) * (Ein.Grid.Step)^2;
Esignal = Esignal*exp(1i*p.Results.phase);

p = real(Esignal);
q = imag(Esignal);

end