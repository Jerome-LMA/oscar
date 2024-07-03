function [Map] = Do_Virtual_Map(G,Power_law,varargin)
%Do_Virtual_map() Create a synthetic map according to a parametrised PSD
%   Grid is the Grid object for the calculations
%   Power_law is a vector of length 8 with the power law, see the function
%   2024 version from M. Lejean and J. Degallaix

p = inputParser;
p.FunctionName = 'Create a virtual map';

% Check if the first argument is a grid
addRequired(p,'G',@(x)isa(x, 'Grid'));

% Check if the PSD law is given
addRequired(p,'Power_law',@(x)isnumeric(x) || ischar(x))

% Check if the phase is given
p.addParameter('phase',[],@(x)isnumeric(x));

% Check if the frequency interval for the PSD_1D generation is given
p.addParameter('bandwidth',[],@(x)isnumeric(x));

% Check if the frequency interval for the PSD_1D generation is given
p.addParameter('porte',false,@(x)islogical(x));

p.parse(G,Power_law,varargin{:})

debug = false;

if ~isnumeric(p.Results.Power_law)
    if strcmp(p.Results.Power_law,'IBF')
        Power_law = [60.13,2.22e+02,2.44,3.24,0.73,2.61]; 
    elseif strcmp(p.Results.Power_law,'Mechanical')
        Power_law = [82.26,3.26e+02,2.49,2.65,1.35,3.54];
    elseif strcmp(p.Results.Power_law,'Coated')
        Power_law = [65.46,3.22e+02,8.90,2.84,0.51,3.70]; 

    else
        Power_law = p.Results.Power_law;
    end
else
    Power_law = p.Results.Power_law;
end

if ~isempty(p.Results.phase)
    Phase = p.Results.phase; % Given Phase
else % Random Phase
    M = zeros(G.Num_point);
    for i = 1 : G.Num_point
        M(i,:) = pi - 2*pi*rand(1,G.Num_point); %generate random number for the row;
        M(:,i) = -M(i,:); % take the same number and put it in the antisymmetric slot
    end
    Phase = M;
end

FFT_Radius_Axis = G.Axis_FFT(G.Half_num_point+1:end); % will take the radius along the diagonal

% ---- Draw the PSD 1D on the new frequency vector ---- %
PSD_1D  = FitFunctionPSD(Power_law,FFT_Radius_Axis);
PSD_1D(1) = PSD_1D(2);

% ---- Reconstruct PSD 2D assuming radial symmetry ---- %
PSD_2D =  interp1(FFT_Radius_Axis,PSD_1D,G.D2_FFT_r,'linear',0);

% ---- Get the map of associated PSD 2D with Phase ---- %
ASD = sqrt(PSD_2D.*G.Length^2/G.Step^2);
Map = ifft2(ifftshift(ASD.*exp(1i*Phase)),'symmetric');
Map = -1*Map; % To discuss 

end