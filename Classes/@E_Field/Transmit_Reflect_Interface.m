function [Eout, varargout] = Transmit_Reflect_Interface(Ein,I1,varargin)
%  Transmit_Reflect_Interface(varargin) Transmit and reflect a E_field through an
%  interface
%  Eout = Transmit_Reflect_Interface(Ein,I1), transmit the E_field Ein
%  across the interface I1
%  [Eout Eref] = Transmit_Reflect_Interface(Ein,I1), transmit the E_field Ein
%  across the interface I1. Eout is the transmitted field and Eref is the reflected field

p  = inputParser;

% Check if the first argument is an E_Field
p.addRequired('Ein', @(x)isa(x, 'E_Field'));

% Check what are we want calculate
p.addRequired('I1', @(x) isa(x, 'Interface') );

p.parse(Ein,I1,varargin{:})

E = Ein;
Inter = I1;

switch nargout
    case 1 % Just compute the transmitted field

        if (E.Refractive_index~=Inter.n1) && (E.Refractive_index~=Inter.n2)
            error('Transmit_Reflect_Interface(): Inconstitrency with the refractive index: the laser light refractive index must match one of the interface refractive index.')
        end

        Eout = E;

        if (E.Refractive_index==Inter.n1)
            PF_Mirror = exp(1i * Eout.k_prop * ((Inter.n2 - Inter.n1)/ Inter.n1) * Inter.surface) .* Inter.mask .* Inter.t;
            Eout.ABCD_q = q_param_transform_ABCD(Inter.ABCD_trans_from_n1,E.ABCD_q);
            Eout = Change_E_n(Eout,Inter.n2);
        else
            % !!! Minus sign on the surface definition
            PF_Mirror = exp(-1i * Eout.k_prop * ((Inter.n1 - Inter.n2)/ Inter.n2) * fliplr(Inter.surface)) .*  fliplr(Inter.mask) .*  Inter.t;
            Eout.ABCD_q = q_param_transform_ABCD(Inter.ABCD_trans_from_n2,E.ABCD_q);
            Eout = Change_E_n(Eout,Inter.n1);
        end
        Eout = Eout .*PF_Mirror;
        varargout = {};

    case 2 % Compute the transmitted AND the reflected field from the interface

        if (E.Refractive_index~=Inter.n1) && (E.Refractive_index~=Inter.n2)
            error('Transmit_E_field(): Inconstitrency with the refractive index: the laser light refractive index must match one of the interface refractive index.')
        end

        Eout = E;
        Eref = E;
        if (E.Refractive_index==Inter.n1)
            PF_Mirror_trans =  exp(1i * Eout.k_prop * ((Inter.n2 - Inter.n1)/ Inter.n1) * Inter.surface)  .* Inter.mask .* Inter.t;
            Eout.ABCD_q = q_param_transform_ABCD(Inter.ABCD_trans_from_n1,E.ABCD_q);
            Eout = Change_E_n(Eout,Inter.n2);

            PF_Mirror_ref = exp(-1i * Eref.k_prop * Inter.surface *2) .* Inter.mask .* Inter.r;
            ABCD_ref = Inter.ABCD_ref_from_n1;
        else
            PF_Mirror_trans =  exp(-1i * Eout.k_prop * ((Inter.n1 - Inter.n2)/ Inter.n2) * fliplr(Inter.surface))  .* fliplr(Inter.mask) .* Inter.t;
            Eout.ABCD_q = q_param_transform_ABCD(Inter.ABCD_trans_from_n2,E.ABCD_q);
            Eout = Change_E_n(Eout,Inter.n1);

            PF_Mirror_ref = exp(1i * Eref.k_prop * fliplr(Inter.surface) *2)  .* fliplr(Inter.mask) .* Inter.r;
            ABCD_ref = Inter.ABCD_ref_from_n2;
        end

        Eout = Eout .* PF_Mirror_trans;

        Eref = Eref .* PF_Mirror_ref;
        Eref.Field = fliplr(Eref.Field);

        Eref.ABCD_q = q_param_transform_ABCD(ABCD_ref,E.ABCD_q);


        if E.Nb_Pair_SB  % if sidebands are present
            for ii=1:E.Nb_Pair_SB
                Eref.SB(ii).Field_lower = fliplr(Eref.SB(ii).Field_lower);
                Eref.SB(ii).Field_upper = fliplr(Eref.SB(ii).Field_upper);
            end
        end

        varargout = {Eref};
end

end


