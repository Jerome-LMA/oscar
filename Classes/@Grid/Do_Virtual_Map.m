function [map] = Do_Virtual_Map(G,Power_law,varargin)
%Do_Virtual_map() Create a synthetic map according to a parametrised PSD
%   Grid is the Grid object for the calculations
%   Power_law is a vector of length 8 with the power law, see the function
%   2024 version from M. Le Jean and J. Degallaix

p = inputParser;
p.FunctionName = 'Create a virtual map';

% Check if the first argument is a grid
addRequired(p, 'G', @(x)isa(x, 'Grid'));

% Check if the PSD law is given
addRequired(p,'Power_law', @(x)isnumeric(x) || ischar(x))

% Check if the phase is given
p.addParameter('phase',[], @(x)isnumeric(x));

p.parse(G,Power_law, varargin{:})

debug = false;

if ~isnumeric(Power_law)
    if strcmp(Power_law, 'IBF')
        Power_law_num = [1.100, 4.728, 43.542, 3.391e+02, 1.636, 3.333, 0.622, 4.016];
    elseif strcmp(Power_law, 'Standard')
        Power_law_num = [11.025, 28.826, 1.445e+02, 2.121, 1.066, 4.641];
    elseif strcmp(Power_law, 'Coated')
        Power_law_num = [99.833, 7.890, 24.510, 3.254e+02, 1.776, 7.165, 1.476, 5.807];
    else
        error('Do_Virtual_Map(): name of the PSD not accepted')
    end
else
    Power_law_num = Power_law;
end

if ~isempty(p.Results.phase)
    Phase = p.Results.phase; % Given Phase
else % Random Phase
    Phase = zeros(G.Num_point);
    for i = 1 : G.Num_point
        Phase(i,:) = pi - 2*pi*rand(1, G.Num_point); % generate random number for the row;
        Phase(:,i) = -Phase(i,:); % put the same numbers in the antisymmetric slots
    end
end

frequency_axis = G.Axis_FFT(G.Half_num_point+1:end); % will take the radius along the diagonal
% Avoid frequency == 0

if(frequency_axis(1,1)==0)
    frequency_axis = frequency_axis(2:end);
end

% ---- Draw the PSD 1D on the new frequency vector ---- %
PSD_1D  = FitFunctionPSD(frequency_axis, Power_law_num);

% ---- Reconstruct PSD 2D assuming radial symmetry ---- %
PSD_2D =  interp1(frequency_axis, PSD_1D, G.D2_FFT_r,'spline',0);

% ---- Get the map of associated PSD 2D with Phase ---- %
ASD = sqrt(PSD_2D / G.Step_sq) * G.Num_point;

%ASD = sqrt(PSD_2D*G.Length^2)/G.Step^2;
map = ifft2(ifftshift(ASD.*exp(1i*Phase)),'symmetric');
map = -1*map; 

end