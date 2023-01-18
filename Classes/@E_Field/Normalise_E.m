function Eout = Normalise_E(Ein, varargin)
% Normalise a E_field
% Use Normalise_E(E_field)    Normalises the E_field to 1W
%        Normalise_E(E_field,P)    Set the E_field to have a
%        power of P Watt
p  = inputParser;

% Check if the first argument is an E_Field
p.addRequired('Ein', @(x)isa(x, 'E_Field'));

% Check if we normalise just the carrier or all
p.addOptional('power', 1, @(x) isnumeric(x)  && (x>=0));

% Check if we normalise just the carrier or all
p.addOptional('include','carrier', @(x)strcmpi(x,'carrier') | strcmpi(x,'all') );

p.parse(Ein,varargin{:})

if p.Results.power == 0
    Eout.Field = complex(zeros(E.Grid.Num_point,Ein.Grid.Num_point));
    if  Ein.Nb_Pair_SB
        for ii=1:Ein.Nb_Pair_SB
            Eout.SB(ii).Field_lower = Eout.Field;
            Eout.SB(ii).Field_upper = Eout.Field;
        end
     end
     return;
end

if strcmp(p.Results.include,'carrier')
    Eout = Ein;
    Eout.Field = Eout.Field / sqrt(Calculate_Power(Ein));
    Eout.Field = Eout.Field * sqrt(p.Results.power);
elseif strcmp(p.Results.include,'all')
    Eout = Ein;
    P_tot = Calculate_Power(Ein, 'include', 'all');
    norm_pow = sqrt(p.Results.power / P_tot);
    
    Eout.Field = Eout.Field * norm_pow;
    
    if Ein.Nb_Pair_SB        
        for ii = 1:Ein.Nb_Pair_SB
            Eout.SB(ii).Field_lower = Eout.SB(ii).Field_lower * norm_pow;
            Eout.SB(ii).Field_upper = Eout.SB(ii).Field_upper * norm_pow;
        end
    end
end      
        
        
        


