function I_Plot(Iin,varargin)
% I_plot(Interface) Plot the surface of an Interface object in an existing
% window

p = inputParser;
p.FunctionName = 'Display an interface object';

% Check if the first argument is an interface
p.addRequired('Iin', @(x)isa(x, 'Interface'));

% Check if the resolution of the grid if given
p.addParameter('diam',[],@(x)isnumeric(x) && x>0);

% Check if the resolution of the grid if given
p.addParameter('zoom',[],@(x)isnumeric(x) && x>0);

% Check if the resolution of the grid if given
p.addParameter('scale',[],@(x)isnumeric(x) && x>0);

p.parse(Iin,varargin{:})

Title_fig = ['Surface profile: ' strrep(inputname(1),'_','\_')];

if ~isempty(p.Results.scale)
    Iin.surface = Iin.surface * p.Results.scale;
end

if isempty(p.Results.diam)
    
    imagesc(Iin.Grid.Axis,Iin.Grid.Axis(end:-1:1),Iin.surface)
    shading interp
    title(Title_fig)
    axis tight
    axis square
    axis xy
    view([0 90])
    colorbar
    
    if ~isempty(p.Results.zoom)
        Axe_limit = Iin.Grid.Length / (2 * p.Results.zoom);
        xlim([-Axe_limit Axe_limit])
        ylim([-Axe_limit Axe_limit])
    end
    
end

if ~isempty(p.Results.diam)
    
    map_ind_cent = Iin.Grid.D2_r > p.Results.diam/2;
    Iin.surface(map_ind_cent) = NaN;
    
    h = imagesc(Iin.Grid.Axis,Iin.Grid.Axis(end:-1:1),Iin.surface);
    shading interp
    title(Title_fig)
    axis tight
    axis square
    axis xy
    view([0 90])
    colorbar
    
    xlim([-p.Results.diam/2 p.Results.diam/2])
    ylim([-p.Results.diam/2 p.Results.diam/2])
    set(h,'AlphaData',~isnan(Iin.surface))
    
end

end