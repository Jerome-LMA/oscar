function [Eout, varargout] = Transmit_Reflect_Interface(varargin)
%  Transmit_Reflect_Interface(varargin) Transmit and reflect a E_field through an
%  interface
%  Eout = Transmit_Reflect_Interface(Ein,I1), transmit the E_field Ein
%  across the interface I1
%  [Eout Eref] = Transmit_Reflect_Interface(Ein,I1), transmit the E_field Ein
%  across the interface I1. Eout is the transmitted field and Eref is the reflected field 

switch nargin
    case {0,1}
        error('Transmit_Reflect_Interface(): Not enough arguments, at least one object E_field and one oject Interface must be given')
        
    case 2
        switch nargout
            case 1 % Just compute the transmitted field
                E = varargin{1};
                Inter = varargin{2};
                
                if (E.Refractive_index~=Inter.n1) && (E.Refractive_index~=Inter.n2)
                    error('Transmit_Reflect_Interface(): Inconstitrency with the refractive index: the laser light refractive index must match one of the interface refractive index.')
                end
                
                
                if isa(E,'E_Field') && isa(Inter, 'Interface')
                    Eout = E;
                    
                    if (E.Refractive_index==Inter.n1)
                        PF_Mirror = exp(1i * Eout.k_prop * ((Inter.n2 - Inter.n1)/ Inter.n1) * Inter.surface) .* Inter.mask .* Inter.t;
                        Eout = Change_E_n(Eout,Inter.n2);
                    else
                        % !!! Minus sign on the surface definition
                        PF_Mirror = exp(-1i * Eout.k_prop * ((Inter.n1 - Inter.n2)/ Inter.n2) * fliplr(Inter.surface)) .*  fliplr(Inter.mask) .*  Inter.t;
                        Eout = Change_E_n(Eout,Inter.n1);
                    end
                    
                    Eout.Field = Eout.Field .*PF_Mirror;
                    
                    if ~isempty(E.Field_SBl) % if sidebands are present
                        Eout.Field_SBl = Eout.Field_SBl  .* PF_Mirror;
                        Eout.Field_SBu = Eout.Field_SBu .* PF_Mirror;
                    end
                    
                    varargout = {};
                    return
                else
                    error('Transmit_Reflect_Interface(): The first argmument must be an object E_Field and the second an interface')
                end
                
            case 2 % Compute the transmitted AND the reflected field from the interface
                E = varargin{1};
                Inter = varargin{2};
                
                if (E.Refractive_index~=Inter.n1) && (E.Refractive_index~=Inter.n2)
                    error('Transmit_E_field(): Inconstitrency with the refractive index: the laser light refractive index must match one of the interface refractive index.')
                end
                
                if isa(E,'E_Field') && isa(Inter, 'Interface')
                    Eout = E;
                    Eref = E;
                    if (E.Refractive_index==Inter.n1)
                        PF_Mirror_trans =  exp(1i * Eout.k_prop * ((Inter.n2 - Inter.n1)/ Inter.n1) * Inter.surface)  .* Inter.mask .* Inter.t;
                        Eout = Change_E_n(Eout,Inter.n2);
                        
                        PF_Mirror_ref = exp(-1i * Eref.k_prop * Inter.surface *2) .* Inter.mask .* Inter.r;
                    else
                        PF_Mirror_trans =  exp(-1i * Eout.k_prop * ((Inter.n1 - Inter.n2)/ Inter.n2) * fliplr(Inter.surface))  .* fliplr(Inter.mask) .* Inter.t;
                        Eout = Change_E_n(Eout,Inter.n1);
                        
                        PF_Mirror_ref = exp(1i * Eref.k_prop * fliplr(Inter.surface) *2)  .* fliplr(Inter.mask) .* Inter.r;
                    end
                    
                    Eout.Field = Eout.Field .* PF_Mirror_trans;
                    
                    Eref.Field = Eref.Field .* PF_Mirror_ref;
                    Eref.Field = fliplr(Eref.Field);
                    
                    if ~isempty(E.Field_SBl) % if sidebands are present
                        Eout.Field_SBl = Eout.Field_SBl  .* PF_Mirror_trans;
                        Eout.Field_SBu = Eout.Field_SBu .* PF_Mirror_trans;
                        
                        Eref.Field_SBl = Eref.Field_SBl .* PF_Mirror_ref;
                        Eref.Field_SBl = fliplr(Eref.Field_SBl);
                        
                        Eref.Field_SBu = Eref.Field_SBu .* PF_Mirror_ref;
                        Eref.Field_SBu = fliplr(Eref.Field_SBu);                        
                    end
                    
                    varargout = {Eref};
                    return
                else
                    error('Transmit_Reflect_Interface(): The first argmument must be an object E_Field and the second an interface')
                end
                
            otherwise
                error('Transmit_Reflect_Interface(): Wrong number of output argument')
        end
    otherwise
        error('Transmit_Reflect_Interface(): invalid number of input arguments, no power calculation is made')
end
end
