function varargout = Weighted_RMS(Iin,varargin)
% Weighted_RMS(Iin,diam) :Calculate the RMS of the Interface Iin weighted by the intensity of the
% beam Ein (an instance of the class E_Field).
% The RMS can also be calculated over a diameter with for example:
% Weighted_RMS(Iin,'diam',0.2) calculate the RMS over a diameter of 0.2 m
% Weighted_RMS(Iin,'E', E_Field(G1,'w0',0.02))

p  = inputParser;

% Check if the first argument is an interface
p.addRequired('Iin', @(x)isa(x, 'Interface'));

% Another parameter could be the diameter
p.addParameter('diam',[],@(x)isnumeric(x) && x>0);

% or weight the fit with the power of a Gaussian beam
p.addParameter('E',[], @(x)isa(x, 'E_Field'));


p.parse(Iin,varargin{:})

Iin = p.Results.Iin;

if ~isempty(p.Results.diam)
    Diam_RMS = p.Results.diam;
    
    map_ind_cent = Iin.Grid.D2_r > Diam_RMS/2;
    Iin.surface(map_ind_cent) = NaN;
    
    rms_w = std(Iin.surface(~isnan(Iin.surface)));

elseif ~isempty(p.Results.E)
    
    Ein = p.Results.E;
    
    Surface = Iin.surface .* Iin.mask;
    
    Beam_intensity = abs(Ein.Field).^2;
    Weighted_avg = sum(sum(Surface.*Beam_intensity)) /  sum(sum(Beam_intensity));
    
    rms_w = sum(sum(Beam_intensity.*(Surface - Weighted_avg).^2)) / sum(sum(Beam_intensity));
    rms_w = sqrt(rms_w);
    
else
    error('Weighted_RMS(): something went wrong, check the second argument')
end

switch nargout
    case 0
        fprintf('Calculated RMS [nm]: %g \n', rms_w*1E9);
    case 1
        varargout{1} = rms_w;
    otherwise
        error(' Weighted_RMS(): Too many output argument')
end

end

