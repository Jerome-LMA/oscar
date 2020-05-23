function [p, q] = Demodulate_SB(varargin)
% [p q] = Demodulate_SB(E1), demodulate the carrier with the sidebdands
% [p q] are the signal in phase and in quadrature

switch nargin
    case 0
        error('Demodulate_SB(): Not enough arguments, at least an object E_field must be given')
    
    case 1
        E = varargin{1};

        if isempty(E.Field_SBl)
           error('Demodulate_SB(): No sidebands fields are present')
        end
       
        
        Esignal = sum( E.Field .* conj(E.Field_SBu) + conj(E.Field) .* E.Field_SBl, 'all') * (E.Grid.Step)^2;
        
        p = real(Esignal);
        q =imag(Esignal);
       
    otherwise
        error('Demodulate_SB(): invalid number of input argument')
end