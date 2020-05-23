function E_plot_SB(varargin)
% Display the 2D amplitude of the sidebands of a E_Field in a new window
% E_plot(E_Field) display the amplitude of the 2 sidebands in a new window
% E_plot(E_Field,zoom) display the amplitudeof the 2 sidebands a new window
% with a factor zoom centered in the middle of the grid

switch nargin
    case 0
        error('Not enough arguments, at least an object E_field must be given')
        return
    case 1
        zoom_plot = 1;
    case 2
        zoom_plot = varargin{2};
    otherwise
        error('Invalid number of input arguments, no plot is made')
        return
end

if isempty(varargin{1}.Field_SBl)
    error('E_plot_SB(): no sidebands are present')
end


limit = varargin{1}.Grid.Length/(2*zoom_plot);

title_var = strrep(inputname(1),'_', '\_');

subplot(1,2,1)
imagesc(varargin{1}.Grid.Axis,varargin{1}.Grid.Axis,abs(varargin{1}.Field_SBl))
shading interp
axis tight
axis square
axis xy
view([0 90])
axis([-limit limit -limit limit])
title(['Lower sideband ' title_var])

subplot(1,2,2)
imagesc(varargin{1}.Grid.Axis,varargin{1}.Grid.Axis,abs(varargin{1}.Field_SBu))
shading interp
axis tight
axis square
axis xy
view([0 90])
axis([-limit limit -limit limit])
title(['Upper sideband ' title_var])

end