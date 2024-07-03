function [PSD_1D,xFreq] = Plot_PSD(I_in,varargin)
%Plot_PSD Plot the 1D PSD of a surface
% Function used to plot (and return) the 1D PSD derived from an object
% Interface
% Adapted from a set of functions made by M. Galimberti in the years
% 2010-2011 at LMA.
% 2023 version from J. Degallaix and M. Lejean

p = inputParser;

% Check if the first argument is an Interface
p.addRequired('I_in', @(x)isa(x, 'Interface'));

% Check if the diameter for the calculation is given
p.addParameter('diam',[],@(x)isnumeric(x) && x>0);

% Check if the 1D PSD is calculated after the sum over one dimension (by default it is radial)
p.addParameter('rect_1D',[],@(x)islogical(x));

% Check if the 1D PSD should be plotted or not
p.addParameter('display',true,@(x)islogical(x));

% Check if it is needed to weight the PSD with a Gaussian beam
p.addParameter('window_Gaussian',[],@(x)isnumeric(x) && x>0);

% Check RMS at the different steps
p.addParameter('RMS_display',false,@(x)islogical(x));

p.parse(I_in,varargin{:})

if isempty(p.Results.rect_1D)
    Calculation_radial = true;
else
    Calculation_radial = ~p.Results.rect_1D;
end
% -------------------------------------------------------------------------------
%            First prepare the map for the calculation

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
map2 = map(ind_cut,ind_cut);
map_grid_r2 = I_in.Grid.D2_square(ind_cut,ind_cut);
N_point = length(ind_cut);

if p.Results.RMS_display
    fprintf('RMS^2 of the cropped surface map: \t %g \n',std(map2(:))^2) %RMS
end

if isempty(p.Results.window_Gaussian)
    % Apply the window and normalize to keep the power constant
    if exist('hann')
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
    m = exp(- 2*map_grid_r2/rad_gauss^2);
end

% radially symmetric Hann window taken from:
% Quantitative characterization of surface
% topography using spectral analysis
% equation 15

cut_axis = I_in.Grid.Axis(ind_cut);
Hann_rad_sym = zeros(length(cut_axis));

for ii = 1:length(cut_axis)
    for jj = 1:length(cut_axis)
        if (cut_axis(ii)^2 + cut_axis(jj)^2) < (cut_axis(end))^2
            Hann_rad_sym(ii,jj) = 1 + cos( (2*pi*sqrt(cut_axis(ii)^2 + cut_axis(jj)^2)) /(2*cut_axis(end)) );
        end
    end
end

m =  Hann_rad_sym;

P_window =  (sum(abs(m).^2, 'all')/N_point^2); % Power of the window

map2 = map2.*m;
map2 = map2 / sqrt(P_window);

if p.Results.RMS_display
    fprintf('RMS^2 of the cropped windowed surface map: \t %g \n',std(map2(:))^2) %RMS
end

%figure(2); imagesc(log(abs(fftshift(fft2(Hann_rad_sym))))); axis square

%------------------------------------------------------------------------
% Calculate the PSD 2D then put it in 1D

% Now can do the PSD 2D
ASD_2D = fftshift(fft2((map2)))*(I_in.Grid.Step^2);  %
PSD_2D = (abs(ASD_2D).^2)/(diam_PSD^2);

if p.Results.RMS_display
    fprintf('RMS^2 of the 2D PSD: \t %g  \n',sum(PSD_2D(:))/(diam_PSD^2)) %RMS
end

xFreq = (0:N_point/2)/(0.5*diam_PSD);
df = max(diff(xFreq));

if ~Calculation_radial % Pass in 1D by integrating along one direction
    PSD_1D_2s = sum(PSD_2D,2)*(df*N_point);
    PSD_1D = (PSD_1D_2s(1:N_point/2+1)+PSD_1D_2s(end:-1:N_point/2));
    xFreq(1) = [];
    PSD_1D = PSD_1D(2:end);
    xFreq = xFreq(1:end)';
    
    if p.Results.RMS_display
        fprintf('RMS^2 of the integrated rectangular 1D PSD: \t %g  \n',2*sum(PSD_1D)/(N_point)) %RMS
    end
    
    label_y_axis = 'Power Spectral Density [m^2/m^-1]';
    
    
else % Pass in 1D by integrating along a constant frequency
    %---- From PSD_2D to Spectral map ----%
    Step = I_in.Grid.Step;
    GridSize = size(PSD_2D,1);
    Vector = 1:1:GridSize;
    
    %---- Radial Integration in frequency Domain ----%
    FFT_Axis =  -1/(2*Step) + (Vector-1)/diam_PSD;
    [FFT_X,FFT_Y] = meshgrid(FFT_Axis);
    FFT_Metrix = sqrt(FFT_X.^2 + FFT_Y.^2);
    FFT_Radius_Axis = FFT_Axis((GridSize/2+1):end);
    dq = mean(diff(FFT_Radius_Axis));
    PSD_1D  = zeros(1,(GridSize/2));
    
    for ii = 1:length(FFT_Radius_Axis)
        ind_ring = find((FFT_Metrix>FFT_Radius_Axis(ii)) & (FFT_Metrix<=FFT_Radius_Axis(ii)+dq));
        PSD_1D(ii) = sum(PSD_2D(ind_ring))/(length(ind_ring)*(Step^2));
    end
    
    if p.Results.RMS_display
        fprintf('RMS^2 of the integrated radial 1D PSD: \t %g  \n',sum(PSD_1D.*FFT_Radius_Axis)/sum(FFT_Radius_Axis(1:end))) %RMS
    end
    
    xFreq = FFT_Radius_Axis(2:end)';
    PSD_1D = PSD_1D(2:end)';
    
    label_y_axis = 'Power Spectral Density [m^3/m^-1]';
    
end

if p.Results.display
    figure(1)
    loglog(xFreq,PSD_1D,'linewidth',2); set(gca,'fontsize',14);grid on
    xlabel('Spatial frequency [1/m]','FontSize',14)
    ylabel(label_y_axis,'FontSize',14)
    axis tight
end


end

