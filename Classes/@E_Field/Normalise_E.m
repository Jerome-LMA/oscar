function Eout = Normalise_E(varargin)
% Normalise a E_field
% Use Normalise_E(E_field)    Normalises the E_field to 1W
%        Normalise_E(E_field,P)    Set the E_field to have a
%        power of P Watt
switch nargin
    case 0
        disp('Normalise_E(): Not enough arguments, at least an object E_field must be given')
        return
    case 1
        E = varargin{1};
        Eout = E;
        Power_amp_inv = 1/sqrt(Calculate_Power(E));
        
        if isgpuarray(E.Field)
            Eout.Field = arrayfun(@times,Eout.Field,Power_amp_inv);
        else
            Eout.Field = Eout.Field * Power_amp_inv;
        end
    case 2
        E = varargin{1};
        Eout = E;
        Eout.Field = Eout.Field / sqrt(Calculate_Power(varargin{1}));
        Eout.Field = Eout.Field * sqrt(varargin{2});
        
        if varargin{2} == 0
            Eout.Field = complex(zeros(E.Grid.Num_point,E.Grid.Num_point));
            if  E.Nb_Pair_SB
                for ii=1:E.Nb_Pair_SB
                    Eout.SB(ii).Field_lower = Eout.Field;
                    Eout.SB(ii).Field_upper = Eout.Field;
                end
                
            end
        end
        
    otherwise
        error('Normalise_E(): invalid number of input arguments, no power calculation is made')
end

