function Eout = Add_sidebands(varargin)
% Eout = Add_sidebands(Ein,mod_freq,mod_index)
% Add a pair of sidebands of frequency: mod_freq with modulating index: mod_index
% This function simulates a phase modulator

switch nargin
    case {0,1,2}
        error('Add_sidebands(): Not enough arguments, at least an object E_field, a frequency and a frequency depth must be given must be given')
        
    case 3
        E = varargin{1};
        mod_freq = varargin{2};
        mod_index = varargin{3};
    
         Eout = E;
         Eout.Frequency_Offset = mod_freq;
    
         % Change the amplitude for the carrier
         Eout.Field = Eout.Field * besselj(0,mod_index);
         
         % add the 2 sidebands field
         Eout.Field_SBl = - Eout.Field * besselj(1,mod_index);
         Eout.Field_SBu = Eout.Field * besselj(1,mod_index);
       
    otherwise
        error('Add_sidebands(): invalid number of input arguments, no sidebands have been created')
        
end