function [Etrans, varargout] = Transmit_Reflect_Optic(varargin)
% Transmit_Reflect_Optic(Ein,Opt,surface)
%  [Eout Eref] =  Transmit_Reflect_Optic((Ein,M1,'AR') or
%  [Eout Eref] =  Transmit_Reflect_Optic((Ein,I1)
% used to transmit and reflect an E_Field, on either an Interface object or
% a mirror Object
% This function overload Transmit_Reflect_Interface() and
% Transmit_Reflect_Mirror()
switch nargin
    case {0,1}
        error('Transmit_Reflect_Optic(): Not enough arguments, at least one object E_field, one interface or mirror must be given')
        
    case 2
        
        if isa(varargin{2},'Interface')
            switch nargout
                case 1 % Just compute the transmitted field
                    Etrans  = Transmit_Reflect_Interface(varargin{1},varargin{2});
                    return
                    
                case 2 % Compute the transmitted AND the reflected field from the interface
                    [Etrans, Eref]  = Transmit_Reflect_Interface(varargin{1},varargin{2});
                    varargout = {Eref};
                    return
            end
            
        elseif isa(varargin{2},'Mirror')
            switch nargout
                case 1 % Just compute the transmitted field
                    Etrans  = Transmit_Reflect_Mirror(varargin{1},varargin{2},'HR');
                    return
                    
                case 2 % Compute the transmitted AND the reflected field from the mirror
                    error('Transmit_Reflect_Optic(): To calculate the reflected field from a mirror, the first surface encountered must be given')
            end
            
        else
            error('Transmit_Reflect_Optic(): The second argument must be a mirror or interface')
            
        end
        
    case 3
        
        if isa(varargin{2},'Interface')
            switch nargout
                case 1 % Just compute the transmitted field
                    Etrans  = Transmit_Reflect_Interface(varargin{1},varargin{2});
                    return
                    
                case 2 % Compute the transmitted AND the reflected field from the interface
                    [Etrans, Eref]  = Transmit_Reflect_Interface(varargin{1},varargin{2});
                    varargout = {Eref};
                    return
            end
            
        elseif isa(varargin{2},'Mirror')
            switch nargout
                case 1 % Just compute the transmitted field
                    Etrans  = Transmit_Reflect_Mirror(varargin{1},varargin{2},varargin{3});
                    return
                    
                case 2 % Compute the transmitted AND the reflected field from the mirror
                    [Etrans, Eref]  = Transmit_Reflect_Mirror(varargin{1},varargin{2},varargin{3});
                    varargout = {Eref};
            end
        else
            error('Transmit_Reflect_Optic(): The second argument must be a mirror or interface')
            
        end
    otherwise
        error('Transmit_Reflect_Mirror(): invalid number of input arguments')
end








end
