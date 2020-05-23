function varargout = Expand_Zernike(Iin,varargin)
% Expand_zernike(): expand a surface on the Zernike polynomials up to a certain
% order
%  Iout = Expand_Zernike(Iin,'Z_order',10,'diam',0.2), take the interface
%  Iin and find the development in Zernike polynomials up to the order 10
%  and over a diameter of 0.2 m. Iout is the new interface as re-created by
%  the Zernike polynomials.
%

p  = inputParser;
p.FunctionName = 'Display an interface object';

% Check if the first argument is an interface
p.addRequired('Iin', @(x)isa(x, 'Interface'));

% Check if the resolution of the grid if given
p.addParamValue('Z_order',2,@(x)isnumeric(x) && x>0);

% Check if the resolution of the grid if given
p.addParamValue('diam',[],@(x)isnumeric(x) && x>0);

p.parse(Iin,varargin{:})
Z_order = p.Results.Z_order;
I_in = p.Results.Iin;


% If the diameter is not given, do the calculations on the diameter of the
% clear aperture.
if ~isempty(p.Results.diam)
    diam_Zernike = p.Results.diam;
else
    max_radius_CA = max(max(I_in.Grid.D2_r(logical(I_in.mask))));
    diam_Zernike = max_radius_CA*2;
end

% Cut the map over the diameter of calculation
Map.inds = find(abs(I_in.Grid.Axis) <= (diam_Zernike/2));
Map.data = Iin.surface(Map.inds,Map.inds);

[Map.data2 Map.coeff_Z] = ZernikeDecomposition(Map.data, Z_order);

% figure(2)
% imagesc(Map.data2); axis square

I_out = Iin;
I_out.surface = zeros(size(I_in.surface));

% Create a new aperture with the same diameter as the one use for the
% Zernike calculation
mask_index = I_in.Grid.D2_r < (diam_Zernike/2);
I_out.mask = zeros(size(I_in.surface));
I_out.mask(mask_index) = 1;

I_out.surface(Map.inds,Map.inds) = Map.data2;

switch nargout
    case 0
        I_plot(I_out,'diam',diam_Zernike)
    case 1
        varargout{1} = I_out;   % Surface recreated by the Zernike
    case 2
        varargout{1} = I_out;      
        varargout{2} = Map.coeff_Z; % Coefficient of the Zernike    
    otherwise
        error('Expand_zernike(): Too many output argument')
end


end

