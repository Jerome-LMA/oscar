function Eout = Add_Sidebands(Ein,Mod_freq,Mod_index,varargin)
% Eout = Add_sidebands(Ein,mod_freq,mod_index)
% Add a pair of sidebands of frequency: Mod_freq with modulating index: Mod_index
% This function simulates a phase modulator

p  = inputParser;
p.FunctionName = 'Add sidebands to a field';

% Check if the first argument is an E_Field
p.addRequired('Ein', @(x)isa(x, 'E_Field'));

% Check if the first argument is an E_Field
p.addRequired('Mod_freq', @(x)isnumeric(x) && x>0);

% Check if the first argument is an E_Field
p.addRequired('Mod_index', @(x)isnumeric(x) && x>0);

p.parse(Ein,Mod_freq,Mod_index)

Eout = Ein;
Eout.SB(Eout.Nb_Pair_SB+1).Frequency_Offset = Mod_freq;
Eout.SB(Eout.Nb_Pair_SB+1).Input_Mod_index = Mod_index;

% Change the amplitude for the carrier
Eout.Field = Eout.Field * besselj(0,Mod_index);

% add the 2 sidebands field
Eout.SB(Eout.Nb_Pair_SB+1).Field_lower = - Eout.Field * besselj(1,Mod_index);
Eout.SB(Eout.Nb_Pair_SB+1).Field_upper = Eout.Field * besselj(1,Mod_index);

% We have added one pair of sidebands
Eout.Nb_Pair_SB = Eout.Nb_Pair_SB + 1;

end