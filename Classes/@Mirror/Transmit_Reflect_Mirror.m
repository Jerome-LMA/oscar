function [Eout, varargout] = transmit_reflect_mirror(varargin)
%  transmit_reflect_mirror(varargin) Transmit and reflect a E_field through
%  a thick mirror. Deal with the mirror as a low finesse FP cavity

%  Eout =  transmit_reflect_mirror(Ein,M1,'HR'), transmit the E_field Ein
%  througth the mirror M1, entering by the HR surface
%  [Eout Eref] =  transmit_reflect_mirror((Ein,M1,'AR'), transmit the E_field Ein
%  across the mirror M1. Eout is the transmitted field and Eref is the reflected field

switch nargin
    case {0,1,2}
        error('transmit_reflect_mirror(): Not enough arguments, at least one object E_field, one Mirror and the first surface meet must be given')
        
    case 3
        
        if isa(varargin{1},'E_Field')
            error('transmit_reflect_mirror(): the first argument must be a E_Field')
        end
        
        if isa(varargin{2},'Mirror')
            error('transmit_reflect_mirror(): the second argument must be a Mirror')
        end
        
        if strcmp(varargin{3},'HR') || strcmp(varargin{3},'AR')
            error('transmit_reflect_mirror(): the third argument must either be the string HR or AR')
        end
        
        E1 = varargin{1};
        M = varargin{2};
        Side_str =  varargin{3};
        
        Eout = E1;
        Field_trans = Normalise_E(E1,0);
        
        if strcmp(Side_str,'HR')
            % i phase shift induced only by the transmission on the HR surface
            
            [Field_tmp,Field_ref] =  Transmit_Reflect_Interface(Eout,M.I_HR);
            
            for qq =1:M.RT_inside
                Field_tmp =   Propagate_E(Field_tmp,M.length_substrate); % From HR to AR
                [Field_tmp_trans, Field_tmp_ref] =  Transmit_Reflect_Interface(Field_tmp,M.I_AR);
                Field_trans = Field_trans + 3*1i*Field_tmp_trans;
                
                Field_tmp = Propagate_E(Field_tmp_ref,M.length_substrate);
                [Field_tmp_trans, Field_tmp_ref] =  Transmit_Reflect_Interface(Field_tmp,M.I_HR);
                
                Field_ref = Field_ref + Field_tmp_trans;
                Field_tmp = Field_tmp_ref;
            end
            
        elseif strcmp(Side_str,'AR')
            
            [Field_tmp Field_ref] =  Transmit_Reflect_Interface(Eout,M.I_AR);
            
            for qq =1:M.RT_inside
                Field_tmp =   Propagate_E(Field_tmp,M.length_substrate); % From AR to HR
                [Field_tmp_trans Field_tmp_ref] =  Transmit_Reflect_Interface(Field_tmp,M.I_AR);
                Field_trans = Field_trans + Field_tmp_trans;
                
                Field_tmp = Propagate_E(Field_tmp_ref,M.length_substrate);
                [Field_tmp_trans Field_tmp_ref] =  Transmit_Reflect_Interface(Field_tmp,M.I_HR);
                
                Field_ref = Field_ref +1i*3* Field_tmp_trans;
                Field_tmp = Field_tmp_ref;
            end
            
        else
            error('transmit_reflect_mirror(), the third argument must be the string AR or HR')
        end
        
        
        switch nargout
            case 1 % Just compute the transmitted field
                Eout  = Field_trans;
                
            case 2 % Compute the transmitted AND the reflected field from the interface
                Eout  = Field_trans;
                varargout = {Field_ref};
                return
                
                
            otherwise
                error('transmit_reflect_mirror(): Wrong number of output argument')
        end
    otherwise
        error('transmit_reflect_mirror(): invalid number of input arguments, no power calculation is made')
end
end
