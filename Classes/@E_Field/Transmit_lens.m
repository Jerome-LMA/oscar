function [Eout] = Transmit_lens(varargin)
%Transmit_lens(E_Field,focal length) Propagate the E_field through a lens
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Transmit_lens(E_in,f_l) pass the field E_in through a
%     thin lens of focal length f_l
%     Transmit_lens(E_in,I1) pass the field E_in through an interface
%     I1 with the wavefront distortion given by k_prop * (n2 - n1) *
%     surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargin
    
    case {0,1}
        disp('Transmit_lens(): not enough arguments, at least an object E_field and a focal length or surface must be given')
        return
    case 2
        Ein = varargin{1};
        
        if isreal(varargin{2})
            f_lens = varargin{2};
            Eout = Ein;
            % Defined the wavefront induced by the lens
            WF_lens = ((2*f_lens) - sign(f_lens)*sqrt((2*f_lens)^2 - Ein.Grid.D2_square))*2;
            PF_lens = exp(1i * WF_lens*Ein.k_prop);
            
            Eout.Field = Ein.Field .* PF_lens;
            
            if ~isempty(Ein.Field_SBl) % if sidebands are present
                Eout.Field_SBl = Ein.Field_SBl  .* PF_lens;
                Eout.Field_SBu = Ein.Field_SBu .* PF_lens;
            end
                       
        elseif isa(varargin{2}, 'Interface')
            I = varargin{2};
            % consider the surface as a thin lens, one surface flat and
            % the other curved.
            Eout = Ein;
            PF_lens = exp(1i * Ein.k_prop * (I.n2 - I.n1) * I.surface) .* I.mask;
            
            Eout.Field = Ein.Field .*  PF_lens;
            
            if ~isempty(Ein.Field_SBl) % if sidebands are present
                Eout.Field_SBl = Ein.Field_SBl  .* PF_lens;
                Eout.Field_SBu = Ein.Field_SBu .* PF_lens;
            end
                       
        else
            disp('Transmit_lens(): The second argument must be a length or an interface')
            return
        end
        
    otherwise
        disp('Transmit_lens(): Invalid number of input arguments')
        return
        
        
end

