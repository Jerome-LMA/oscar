function [PSD_1D,freq] = Plot_PSD(I_in,varargin)
%Plot_PSD Plot the 1D PSD of a surface
% Function used to plot (and return) the 1D PSD derived from an object
% Interface
% Adapted from a set of functions made by M. Galimberti in the years
% 2010-2011 at LMA.

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

P_window =  (sum(abs(m).^2, 'all')/N_point^2); % Power of the window

map2 = map2.*m;
map2 = map2 / sqrt(P_window);

%std(map2(:))^2

%------------------------------------------------------------------------
% Calculate the PSD 2D then put it in 1D

% Now can do the PSD 2D
ASD_2D = fft2(ifftshift(map2)) * I_in.Grid.Step^2;  % !! low frequency at the corner now
PSD_2D = abs(ASD_2D).^2 / diam_PSD^2;

% figure(1)
% imagesc(abs(fftshift(fftshift(ASD_2D,1),2))); axis square



Vec_f = (0: N_point/2) * 1/diam_PSD;
df = max(diff(Vec_f));

if ~Calculation_radial % Pass in 1D by integrating along one direction
    PSD_1D_2s = sum(PSD_2D,2) * df;
    PSD_1D = PSD_1D_2s(1:N_point/2 +1) + PSD_1D_2s(end:-1:N_point/2);
    
else % Pass in 1D by integrating along a constant radius
    rad_au = size(PSD_2D,1)/2;
    [rad_auX,rad_auY] = meshgrid(-rad_au: rad_au-1);
    rad_au2D =sqrt(rad_auX.^2 + rad_auY.^2);

    PSD_1D = zeros(1,rad_au);
    PSD_2D = fftshift(PSD_2D); % Bring back the low frequency in the middle of the map
    
    for kk=1:rad_au
        ind_area = rad_au2D >= kk - 1  & rad_au2D < kk;
        PSD_1D(kk) = sum(PSD_2D(ind_area));
    end
    PSD_1D = PSD_1D * df;
    Vec_f(1) = [];
    
end

freq = Vec_f;

if p.Results.display
    loglog(Vec_f,PSD_1D,'linewidth',2); set(gca,'fontsize',14);grid on
    xlabel('Spatial frequency [1/m]','FontSize',14)
    ylabel('Power Spectral Density [m^2/m^-1]','FontSize',14)
    axis tight
end



end

