function E_plot(Ein,varargin)
% Display the 2D amplitude of a E_Field in a new window
% E_plot(E_Field) display the amplitude of the field in a new window
% E_plot(E_Field,zoom) display the amplitude of the field in a new window
% with a factor zoom centered in the middle of the grid

p  = inputParser;
p.FunctionName = 'Display an E_Field object';

% Check if the first argument is an E_Field
p.addRequired('Ein', @(x)isa(x, 'E_Field'));

% Check if we display the amplitude or intensity
p.addOptional('display','amplitude', @(x)strcmpi(x,'amplitude') | ...
    strcmpi(x,'intensity'));

% Check if the scale is linear or logarithmic
p.addOptional('scale','linear', @(x)strcmpi(x,'linear') | ...
    strcmpi(x,'log'));

% Check if we want to display the beam or its divergence
p.addOptional('domain','space', @(x)strcmpi(x,'space') | ...
    strcmpi(x,'angle'));

p.addOptional('no_axis',false, @(x)islogical(x) );

p.parse(Ein,varargin{:})

title_str = strrep(inputname(1),'_', '\_');


if strcmp(p.Results.display,'amplitude')
    if strcmp(p.Results.domain,'space')
        todisplay = abs(Ein.Field); % display the amplitude
        ax_scale = Ein.Grid.Axis;
        title1 = 'Amplitude';
    else
        todisplay = abs(fftshift(fft2(Ein.Field))); % display the angle
        ax_scale_tmp = Ein.Grid.Axis_FFT * Ein.Wavelength;
        ax_scale = ax_scale_tmp(end:-1:1);
        title1 = 'Angular direction';
    end
else
    todisplay = abs(Ein.Field).^2; % display the intensity
    ax_scale = Ein.Grid.Axis;
    title1 = 'Intensity';
end


if strcmp(p.Results.scale,'linear')
    imagesc(ax_scale,ax_scale,todisplay)
    title([title1 ' of the electric field: ' title_str])
else
    imagesc(ax_scale,ax_scale,log(todisplay))
    title([title1 ' of the electric field log scale: ' title_str])
end

shading interp

axis equal
axis tight
view([0 90])

if p.Results.no_axis
    axis off
    title('')
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
end


end


