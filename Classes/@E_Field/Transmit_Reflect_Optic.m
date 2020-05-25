function [Etrans, varargout] = Transmit_Reflect_Optic(Ein,Min,side,varargin)
% Transmit_Reflect_Optic(Ein,Opt,surface)
%  [Eout Eref] =  Transmit_Reflect_Optic((Ein,M1,'AR') or
%  [Eout Eref] =  Transmit_Reflect_Optic((Ein,I1)
% used to transmit and reflect an E_Field, on either an Interface object or
% a mirror Object
% This function overload Transmit_Reflect_Interface() and
% Transmit_Reflect_Mirror()

p  = inputParser;

% Check if the first argument is an E_Field
p.addRequired('Ein', @(x)isa(x, 'E_Field'));

% Check if the first argument is an E_Field
p.addRequired('Min', @(x) (isa(x, 'Mirror') | isa(x, 'Interface')));

% Check what are we want calculate
p.addOptional('side', @(x)strcmpi(x,'HR') | strcmpi(x,'AR') );

p.parse(Ein,Min,varargin{:})

% if the Min is a 'Mirror' and you want to calculate the reflection, the
% third arugment is mandatory

if isa(Min,'Interface')
    switch nargout
        case 1 % Just compute the transmitted field
            Etrans  = Transmit_Reflect_Interface(Ein,Min);
            return
            
        case 2 % Compute the transmitted AND the reflected field from the interface
            [Etrans, Eref]  = Transmit_Reflect_Interface(Ein,Min);
            varargout = {Eref};
            return
    end
    
else
    switch nargout
        case 1 % Just compute the transmitted field
            Etrans  = Transmit_Reflect_Mirror(Ein,Min,'HR');
            return
            
        case 2 % Compute the transmitted AND the reflected field from the mirror
            if ~(exist('side','var') == 1)
                error('Transmit_Reflect_Optic(): To calculate the reflected field from a mirror, the first surface encountered must be given')
            end
            
            [Etrans, Eref]  = Transmit_Reflect_Mirror(Ein,Min,side);
            varargout = {Eref};
    end
    
    
end

end
