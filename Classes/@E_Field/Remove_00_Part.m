function Eout = Remove_00_Part(Ein,varargin)
%Remove_00_part Remove the fundamental Gaussian part from an E_Field
%   Return a E_Field where the fundamental Gaussian part has been
%   substracted. Very useful to highligh the distortion of a beam
%   Eout is normalised by setting the power of the input beam to 1
%  The substraction of 2 beams can also be achieved Remove_00_part(E1,E2)
%  !! ONLY WORK FOR THE CARRIER

p  = inputParser;
p.FunctionName = 'Display an E_Field object';

% Check if the first argument is an E_Field
p.addRequired('Ein', @(x)isa(x, 'E_Field'));

% Check if a second E_field is given 
p.addOptional('E_sub',[], @(x)isa(x, 'E_Field'));

p.parse(Ein,varargin{:})

Ein = p.Results.Ein;
Esub = p.Results.E_sub;

Ein.Nb_Pair_SB = 0;

Ein = Normalise_E(Ein,'Power',1);

% First, fit the beam with a TEM00

if isempty(Esub)
    [E00_w, E00_r] = Fit_TEM00(Ein);
    E00 = E_Field(Ein.Grid,'w',E00_w,'R',E00_r,'Wavelength',Ein.Wavelength);
else
    E00 = Normalise_E(Esub,'Power',1);
end

diff_phase = angle(Calculate_Overlap(Ein,E00));

% To minimize the output, be sure the 2 beams interfer destructively
Eout = Ein - E00 * exp(1i*diff_phase);

E_Plot(Eout)

end

