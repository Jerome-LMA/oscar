function [PSD_1D,xFreq] = Plot_PSD(I_in,varargin)
%Plot_PSD Plot the 1D PSD of a surface
% Function used to plot (and return) the 1D PSD derived from an object
% Interface
% Adapted from a set of functions made by M. Galimberti in the years
% 2010-2011 at LMA.
% 2024 version from J. Degallaix and M. Lejean

p = inputParser;

% Check if the first argument is an Interface
p.addRequired('I_in', @(x)isa(x, 'Interface'));

% Check if the diameter for the calculation is given
p.addParameter('diam',[],@(x)isnumeric(x) && x>0);

% Check if the 1D PSD is calculated after the sum over one dimension (by default it is radial)
p.addParameter('radial',true,@(x)islogical(x));

% Check if the 1D PSD should be plotted or not
p.addParameter('display',true,@(x)islogical(x));

% Check if it is needed to weight the PSD with a Gaussian beam
p.addParameter('window_Gaussian',[],@(x)isnumeric(x) && x>0);

% Check RMS at the different steps
p.addParameter('RMS_display',false,@(x)islogical(x));

p.parse(I_in,varargin{:})
Calculation_radial = p.Results.radial;


% -------------------------------------------------------------------------------
%            Preprocessing the map for the calculation

% Remove the NaN if any
I_in.surface(isnan(I_in.surface)) = 0;

% Calculation only within a diameter or within the diameter given by the mirror
% aperture
if ~isempty(p.Results.diam)
    diam_PSD = p.Results.diam;
else
    max_radius_CA = max(max(I_in.Grid.D2_r(logical(I_in.mask))));
    diam_PSD = max_radius_CA*2;
end

map = I_in.surface;

% Remove the average over the central part
ind_in = find(I_in.Grid.D2_r <= diam_PSD/2);
map = map - mean(map(ind_in));

% Add 0 outside the central part
mask = zeros(size(map));
mask(ind_in) = 1;
map = map .* mask;

% Cut the map around the useful part (a square of length diam_PSD)
ind_cut = find(abs(I_in.Grid.Axis) <= diam_PSD/2);
map = map(ind_cut,ind_cut);
N_point = length(ind_cut);

% Redefine the grid for the calculation and calculate the frequency axis
G0 = Grid(N_point,diam_PSD);
freq_axis = G0.Axis_FFT((N_point/2+1):end); % take positive frequency
df = max(diff(freq_axis));

%figure(102);imagesc(map)

if p.Results.RMS_display
    fprintf('RMS^2 of the cropped surface map: \t %g \n',std(map(:))^2) %RMS !! we have added 0 in the corner
end

if isempty(p.Results.window_Gaussian)
    % Apply the window and normalize to keep the power constant
    if exist('hann','builtin')
        w = hann(N_point);
    else
        w = 0.5 * (1 - cos(2*pi*(0:N_point-1)'/(N_point-1)));
        %w = ones(1,N_point);
    end
    %  w = bohmanwin(N_point);
    %  chebwin
    m = w(:)*w(:).'; % Outer product window
else
    rad_gauss = p.Results.window_Gaussian;
    m = exp(- 2*G0.D2_square/rad_gauss^2);
end

P_window =  (sum(abs(m).^2, 'all')/N_point^2); % Power of the window

map = map.*m;
map = map/sqrt(P_window);


if p.Results.RMS_display
    fprintf('RMS^2 of the cropped windowed surface map: \t %g \n',std(map(:))^2) %RMS
end

%figure(103);imagesc(map)


%------------------------------------------------------------------------
% Calculate the PSD 2D

FFT_2D_norm = fftshift(fft2(map))/N_point;
FFT_2D_norm_sq = abs(FFT_2D_norm).^2;

% According to the Parseval's theorem, the 2 following quantities must be
% equal:
%sum(FFT_2D_norm_sq(:))
%sum(map.^2,'all')

if p.Results.RMS_display
    fprintf('RMS^2 of the 2D FFT squared: \t %g  \n',sum(FFT_2D_norm_sq(:))/N_point^2) %RMS
end

PSD_2D = FFT_2D_norm_sq*I_in.Grid.Step_sq; 
% normalisation to match the equation RMS^2 = int int PSD_2D(fx, fy) * dfx * dfy
% dfx, dfy = 1/L

if p.Results.RMS_display
    fprintf('RMS^2 derived from the 2D PSD: \t %g  \n',sum(PSD_2D(:))*df^2) %RMS
end


%------------------------------------------------------------------------
% Converting the PSD 2D to PSD 1D

if ~Calculation_radial % Pass in 1D by integrating along one direction
    label_y_axis = 'Power Spectral Density [m^2/m^-1]';

    PSD_1D_2s = sum(PSD_2D,2)*df;
    PSD_1D = (PSD_1D_2s(1:N_point/2)+PSD_1D_2s(end:-1:N_point/2+1)); % sum with the positive and negative frequency
     
    PSD_1D = PSD_1D(2:end);
    freq_axis = freq_axis(end:-1:2)';
    
    if p.Results.RMS_display
        fprintf('RMS^2 of the integrated rectangular 1D PSD: \t %g  \n',sum(PSD_1D)*df) %RMS
    end
    

else % Pass in 1D by averaging along rings at constant frequency
    label_y_axis = 'Power Spectral Density [m^3/m^-1]';
    
    freq_coordinates = sqrt(G0.D2_FFT_square);
    
    PSD_1D = zeros(1,G0.Num_point/2);
    for ii = 1:length(freq_axis)
        ring = (freq_coordinates>freq_axis(ii)) & (freq_coordinates<=freq_axis(ii)+df); %
        PSD_1D(ii) = mean(PSD_2D(ring));
    end

    
    %  Dropping the first point because freq = 0 is not physical.
    freq_axis = freq_axis(2:end)';
    PSD_1D = PSD_1D(2:end)';

    if p.Results.RMS_display
        fprintf('RMS^2 of the integrated radial 1D PSD: \t %g  \n' ...
            ,    2*pi*sum(PSD_1D.*freq_axis)*df); 
        % if purely random map, correct for the factor 4/pi to take into
        % account the missing corner when the taking the PSD 1D
    end
end

if p.Results.display
    figure(1)
    loglog(freq_axis,PSD_1D,'linewidth',2); 
    set(gca,'fontsize',14, 'XScale', 'log', 'YScale','log');
    xlabel('Spatial frequency [1/m]','FontSize',14)
    ylabel(label_y_axis,'FontSize',14)
    grid on
    axis tight
end

end

