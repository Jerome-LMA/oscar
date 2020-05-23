function varargout = Weighted_RoC(Iin,Ein)
% Fit_Surface_beam() find the radius of curvature as seens by the laser
% beam (another way to say it: the fit of the radius of curvature is weighted by the beam intensity).
% The tilt is also taken into account during the fit

p  = inputParser;

% Check if the first argument is an interface
p.addRequired('Iin', @(x)isa(x, 'Interface'));

% Check if the second argument is an E_Dield
p.addRequired('Ein', @(x)isa(x, 'E_Field') || (isnumeric(x) && x>0));

p.parse(Iin,Ein)

Iin = p.Results.Iin;
Ein = p.Results.Ein;

if isa(Ein, 'E_Field')   % An E_field is entered for the beam weighting
    
    % Do the fit over the diameter of the mask
    max_radius_CA = max(max(Iin.Grid.D2_r(logical(Iin.mask))));
    
    Map.mask_index = find(Iin.Grid.D2_r < max_radius_CA); % Find the right index
    
    Map.Gridv(:,1) = Iin.Grid.D2_X(Map.mask_index);
    Map.Gridv(:,2) = Iin.Grid.D2_Y(Map.mask_index);
    Map.surface = Iin.surface(Map.mask_index);
    
    Beam_intensity = abs(Ein.Field);
    Map.weight = Beam_intensity(Map.mask_index); % will be squared later in the fir
    
    % Do the fit of the tilt and curvature
    func_curv2 = @(c,xdata) c(1)*(xdata(:,1).^2 +xdata(:,2).^2) + c(2) +...
        c(3)*xdata(:,1) + c(4)*xdata(:,2);
    
    %func_curv2 = @(c,xdata) c(1)*(xdata(:,1).^2 +xdata(:,2).^2) + c(2) 
    
    c0 = [-1/(2*1E3) 0 0 0];
    %c0 = [-1/(2*1E3) 0];

    options = optimset('Display','off','MaxFunEvals',1E6,'TolFun',1E-16,'DiffMinChange',1E-14);
    
    [Map.fit_para2,~,~,~,~] = lsqnonlin(@(p)Map.weight.*(func_curv2(p,Map.Gridv)-Map.surface),c0,[],[],options);
    
    RoC_fitted = 1/(2*Map.fit_para2(1));
    
else % A diameter is entered, so the fit will be do over a diameter
    Diam_RoC = Ein;
    
    % Select the right point
    index_fit = find(Iin.Grid.D2_r < Diam_RoC /2);
    
    Map.Gridv(:,1) = Iin.Grid.D2_X(index_fit);
    Map.Gridv(:,2) = Iin.Grid.D2_Y(index_fit);
    Map.surface = Iin.surface(index_fit);
    
    % Do the fit of the tilt and curvature
    func_curv2 = @(c,xdata) c(1)*(xdata(:,1).^2 +xdata(:,2).^2) + c(2) +...
        c(3)*xdata(:,1) + c(4)*xdata(:,2);
    
    c0 = [-1/(2*1E3) 0 0 0];
    options = optimset('Display','off','MaxFunEvals',1E6,'TolFun',1E-16,'DiffMinChange',1E-14);
    
    [fit_para,~,~,~,output] = lsqcurvefit(func_curv2,c0,Map.Gridv,Map.surface,[],[],options);
    RoC_fitted = 1/(2*fit_para(1));
    
end

switch nargout
    case 0
        fprintf('Calculated RoC [m]: %g \n',RoC_fitted);
    case 1
        varargout{1} = RoC_fitted;
    otherwise
        error('Fit_TEM00(): Too many output argument')
end

end

