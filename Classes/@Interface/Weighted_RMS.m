function varargout = Weighted_RMS(Iin,diam)
% Weighted_RMS(Iin,diam) :Calculate the RMS of the Interface Iin weighted by the intensity of the
% beam Ein (an instance of the class E_Field).
% The RMS can also be calculated over a diameter with for example:
% Weighted_RMS(Iin,0.2) calculate the RMS over a diameter of 0.2 m

p  = inputParser;

% Check if the first argument is an interface
p.addRequired('Iin', @(x)isa(x, 'Interface'));

% Check if the second argument is an E_Dield
p.addRequired('diam', @(x) isa(x, 'E_Field') || (isnumeric(x) && x>0) );

p.parse(Iin,diam)

Iin = p.Results.Iin;

if isa(p.Results.diam, 'E_Field')   % An E_field is entered
    Ein = p.Results.diam;
    
    Surface = Iin.surface .* Iin.mask;
    
    Beam_intensity = abs(Ein.Field).^2;
    Weighted_avg = sum(sum(Surface.*Beam_intensity)) /  sum(sum(Beam_intensity));
    
    rms_w = sum(sum(Beam_intensity.*(Surface - Weighted_avg).^2)) / sum(sum(Beam_intensity));
    rms_w = sqrt(rms_w);
    
else % A diameter is entered
    Diam_RMS = p.Results.diam;
    
    map_ind_cent = Iin.Grid.D2_r > Diam_RMS/2;
    Iin.surface(map_ind_cent) = NaN;

    rms_w = std(Iin.surface(~isnan(Iin.surface)));   
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

