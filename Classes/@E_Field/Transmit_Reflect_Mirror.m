function [Eout, varargout] = Transmit_Reflect_Mirror(varargin)
%  Transmit_Reflect_Mirror(varargin) Transmit and reflect a E_field through
%  a thick mirror. Deal the mirror as a low finesse FP cavity

%  Eout =  Transmit_Reflect_Mirror((Ein,M1,'HR'), transmit the E_field Ein
%  througth the mirror M1, entering by the HR surface
%  [Eout Eref] =  Transmit_Reflect_Mirror((Ein,M1,'AR'), transmit the E_field Ein
%  across the mirror M1. Eout is the transmitted field and Eref is the reflected field

switch nargin
    case {0,1,2}
        error('Transmit_Reflect_Mirror(): Not enough arguments, at least one object E_field, one Mirror and the first surface meet must be given')
        
    case 3
        
        if ~isa(varargin{1},'E_Field')
            error('Transmit_Reflect_Mirror(): the first argument must be a E_Field')
        end
        
        if ~isa(varargin{2},'Mirror')
            error('Transmit_Reflect_Mirror(): the second argument must be a Mirror')
        end
        
        if ~(strcmp(varargin{3},'HR') || strcmp(varargin{3},'AR'))
            error('Transmit_Reflect_Mirror(): the third argument must either be the string HR or AR')
        end
        
        E1 = varargin{1};
        M = varargin{2};
        Side_str =  varargin{3};
        
        Eout = E1;
        Field_trans = Normalise_E(E1,0);
        
        if ~isempty(E1.Field_SBl)
        error('Transmit_Reflect_Mirror(): reflection transmission, not yet implemented for the SB')
        end
        
        if strcmp(Side_str,'HR')
            % i phase shift induced only by the transmission on the HR surface
            
            [Field_tmp, Field_ref] =  Transmit_Reflect_Interface(Eout,M.I_HR);
                        
            for qq =1:M.RT_inside
          
                Field_tmp =   Propagate_E(Field_tmp,M.length_substrate); % From HR to AR
                [Field_tmp_trans, Field_tmp_ref] =  Transmit_Reflect_Interface(Field_tmp,M.I_AR);
                 
                %Calculate_power(Field_tmp_trans)
                %Field_trans = Field_trans + (-1i)*Field_tmp_trans;
                Field_trans = Field_trans + Field_tmp_trans;
         
                
                Field_tmp = Propagate_E(Field_tmp_ref,M.length_substrate);
                [Field_tmp_trans, Field_tmp_ref] =  Transmit_Reflect_Interface(Field_tmp,M.I_HR);
                
                
                
                Field_ref = Field_ref + Field_tmp_trans;
                Field_tmp = Field_tmp_ref;
            end
            
        else
            % Start with the AR surface
            [Field_tmp, Field_ref] =  Transmit_Reflect_Interface(Eout,M.I_AR);
                       
            for qq =1:M.RT_inside
                Field_tmp =   Propagate_E(Field_tmp,M.length_substrate); % From AR to HR

                
                [Field_tmp_trans, Field_tmp_ref] =  Transmit_Reflect_Interface(Field_tmp,M.I_HR);
                %Field_trans = Field_trans + (-1i)*Field_tmp_trans;
                Field_trans = Field_trans + Field_tmp_trans;
                
                Field_tmp = Propagate_E(Field_tmp_ref,M.length_substrate);
                [Field_tmp_trans, Field_tmp_ref] =  Transmit_Reflect_Interface(Field_tmp,M.I_AR);
                
                Field_ref = Field_ref + Field_tmp_trans;               
                Field_tmp = Field_tmp_ref;
            end
           
        end
        
        
        switch nargout
            case 1 % Just compute the transmitted field
                Eout  = Field_trans;
                
            case 2 % Compute the transmitted AND the reflected field from the interface
                Eout  = Field_trans;
                varargout = {Field_ref};
                return
                
                
            otherwise
                error('Transmit_Reflect_Mirror(): Wrong number of output argument')
        end
    otherwise
        error('Transmit_Reflect_Mirror(): invalid number of input arguments, no power calculation is made')
end
end
