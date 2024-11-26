function [Eout] = Transmit_Lens(Ein,f_L,varargin)
%Transmit_lens(E_Field,focal length) Propagate the E_field through a lens
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Transmit_lens(E_in,f_l) pass the field E_in through a
%     thin lens of focal length f_l
%     Transmit_lens(E_in,I1) pass the field E_in through an interface
%     I1 with the wavefront distortion given by k_prop * (n2 - n1) *
%     surface. A positive RoC of I1 will give a divergent lens
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p  = inputParser;

% Check if the first argument is an E_Field
p.addRequired('Ein', @(x)isa(x, 'E_Field'));

% Check what are we want calculate
p.addRequired('f_L', @(x) (isreal(x) || isa(x, 'Interface')));

p.parse(Ein,f_L,varargin{:})

if isreal(f_L)
    f_lens = f_L;
    % Defined the wavefront induced by the lens
    WF_lens = ((2*f_lens) - sign(f_lens)*sqrt((2*f_lens)^2 - Ein.Grid.D2_square))*2;
    PF_lens = exp(1i * WF_lens*Ein.k_prop);
    
    Eout = Ein .* PF_lens;
    Eout.ABCD_q = [1 0;-1/ f_L 1] * Ein.ABCD_q;

    
else
    % consider the surface as a thin lens, one surface flat and
    % the other curved given by the interface.
    
    Eout = Transmit_Reflect_Interface(Ein,f_L);  
    Eout = Change_E_n(Eout,Ein.Refractive_index);
    
end


end

