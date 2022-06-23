function Eout = Add_Sidebands(Ein,varargin)
% Eout = Add_sidebands(Ein,mod_freq,mod_index)
% Add a pair of sidebands of frequency: Mod_freq with modulating index: Mod_index
% This function simulates a phase modulator

p  = inputParser;
p.FunctionName = 'Add sidebands to a field';

% Check if the first argument is an E_Field
p.addRequired('Ein', @(x)isa(x, 'E_Field'));

% Input the modulation frequency in Hertz
p.addParameter('Mod_freq',1E6,@(x)isnumeric(x) && x>0);

% Input the modulation depth for the first order
p.addParameter('Mod_depth',0.1 ,@(x)isnumeric(x) && x>0);

% Input the maximum order to take into account (i.e. how many harmonics),
% only one SB pair is order 1.
p.addParameter('Mod_order',1 ,@(x)isnumeric(x) && x>0);


p.parse(Ein,varargin{:})

Max_mod_order = round(p.Results.Mod_order);
Mod_index = p.Results.Mod_depth;

Eout = Ein;
Eout.SB(Eout.Nb_Pair_SB+1).Frequency_Offset = p.Results.Mod_freq;
Eout.SB(Eout.Nb_Pair_SB+1).Input_Mod_index = Mod_index;

% Change the amplitude for the carrier
Eout.Field = Eout.Field * besselj(0,Mod_index);

% add the 2 sidebands field
for i=1:Max_mod_order
    % We are adding one pair of sidebands
    Eout.Nb_Pair_SB = Eout.Nb_Pair_SB + 1;
    Eout.SB(Eout.Nb_Pair_SB).Field_lower = - Eout.Field * besselj(i,Mod_index);
    Eout.SB(Eout.Nb_Pair_SB).Field_upper = Eout.Field * besselj(i,Mod_index);
end

end
