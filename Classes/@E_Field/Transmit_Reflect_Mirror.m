function [Eout, varargout] = Transmit_Reflect_Mirror(Ein,Min,side,varargin)
%  Transmit_Reflect_Mirror(varargin) Transmit and reflect a E_field through
%  a thick mirror. Deal the mirror as a low finesse FP cavity

%  Eout =  Transmit_Reflect_Mirror((Ein,M1,'HR'), transmit the E_field Ein
%  througth the mirror M1, entering by the HR surface
%  [Eout Eref] =  Transmit_Reflect_Mirror((Ein,M1,'AR'), transmit the E_field Ein
%  across the mirror M1. Eout is the transmitted field and Eref is the reflected field

p  = inputParser;

% Check if the first argument is an E_Field
p.addRequired('Ein', @(x)isa(x, 'E_Field'));

% Check if the first argument is an E_Field
p.addRequired('Min', @(x)isa(x, 'Mirror'));

% Check what are we want calculate
p.addRequired('side', @(x)strcmpi(x,'HR') | strcmpi(x,'AR') );

p.parse(Ein,Min,side,varargin{:})


E1 = Ein;
M = Min;
Side_str = side;

Eout = E1;
Field_trans = Normalise_E(E1,0);



if strcmp(Side_str,'HR')
    % i phase shift induced only by the transmission on the HR surface
    
    [Field_tmp, Field_ref] =  Transmit_Reflect_Interface(Eout,M.I_HR);
    
    for qq =1:M.RT_inside
        
        Field_tmp = Propagate_E(Field_tmp,M.length_substrate); % From HR to AR
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

end
