function Iout = Add_Map(Iin,map_loaded,varargin)
%     Add_map() Load a map and add it to a surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Iout = Add_map(Iin,filename,'reso',1E-3) add the map found in the file
%    given by 'filename' and add it to the interface Iin. 'reso' is
%    the resolution of one pixel of the loaded map. The map must be a
%    square matrix. filename is a string if you want to load a map or a
%    variable name is the map is already in the workspace
%
%    If the map to add has a cylindrical symmetry a 2 columns table could
%    be enough. The resolution is then not necessary, the first column is
%    the radius and the second column is the surface height.
%
%     Iout = Add_map(Iin,filename,'reso',1E-3, 'scale', 2) 'scale' add a
%     scaling to the map, so in that case the map is multiply by 2
%
%     Iout = Add_map(Iin,filename,'reso',1E-3, 'RMS', 5E-9) 'RMS' scale the
%     RMS of map to 5 nm
%
%     Iout = Add_map(Iin,filename,'reso',1E-3, 'rotate', 1) 'rotate' is
%     used if we want the map to rotated by 90 degree. 'rotate',2 to
%     rotate the map by 180 degree
%
%     Iout = Add_map(Iin,filename,'reso',1E-3, 'rotate', 1,'remove_tilt',0.150)
%     Remove the tilt / piston from the map. The argument after is the
%     diameter over which it is removed.
%
%     Iout = Add_map(Iin,filename,'reso',1E-3,'remove_tilt',0.150,'shift',[0.05 0.01])
%     Center the map at the coordinates given (in m). The tilt/focus is
%     removed on this new center. Rotation will happen after this
%     centering.
%
%     Iout = Add_map(Iin,filename,'reso',1E-3, 'rotate', 1,'remove_tilt_focus',0.150)
%     Same as above but also remove the curvature term
%
%
%     Positive value in the map brings the surface closer to the incoming beam
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%
p  = inputParser;
p.FunctionName = 'Load a map';

% Check if the first argument is an interface
p.addRequired('Iin', @(x)isa(x, 'Interface'));

% Check if the second argument is a file to load or a matrix
p.addRequired('map_loaded', @(x)(ischar(x) || ismatrix(x)));

% Check if the resolution of the grid if given
p.addParameter('reso',[],@(x)isnumeric(x) && x>0);

% Check if there is a scaling factor for the map
p.addParameter('scale',[],@(x)isnumeric(x));

% Check if the RMS value of the map must be fixed
p.addParameter('RMS',[],@(x)isnumeric(x) && x>0);

% Check if we should rotate the map by 90 deg and how many time
p.addParameter('rotate',0,@(x)isnumeric(x));

% Check if we should remove the tilt over a certain diameter (in meter)
p.addParameter('remove_tilt',[],@(x)isnumeric(x) && x>0);

% Check if we should remove the tilt and curvarture over a certain diameter (in meter)
p.addParameter('remove_tilt_focus',[],@(x)isnumeric(x) && x>0);

% Check if we should shift the center of the map (in meter)
p.addParameter('shift',[0 0],@(x)isnumeric(x));

% Display or not the results of the fit in the command line
p.addParameter('verbose',true,@(x)isa(x,'logical'));

p.parse(Iin,map_loaded,varargin{:})

%p.Results

Iout = Iin;

if ischar(map_loaded)   % Check if the argument is a filename
    fprintf('Loading the map from file: %s \n',map_loaded)
    map.loaded = load(map_loaded);
elseif ismatrix(map_loaded)   % Check if the argument is a variable
    map.loaded =map_loaded;
else
    error('Add_map(): the second argument must be a filename or a variable')
end

% Check if the matrix is square
[m,n] = size(map.loaded);

if n~=2 % so we do not have 2 columns, but we have a matrix
    if (m~=n) % the map is not square, extract a cutted
        map.loaded = Make_square(map.loaded);
        disp('Add_Map(): the loaded map is made square')
    end
end

if (m==n)     % The matrix is square
    
    % Rescale
    if  ~isempty(p.Results.scale)
        map.loaded = map.loaded * p.Results.scale;
    end
    
    if isempty(p.Results.reso)
        error('Add_map(): Since the map is square, the parameters reso must be given')
    else
        map.res = p.Results.reso;
    end
    
    map.nb_point = m;
    % Create the grid for the loaded map
    map.grid_size = map.nb_point*map.res;
    map.Grid_axis = -map.grid_size/2 + map.res/2 + ((1:map.nb_point)-1)*map.res;
    [map.Grid_X,map.Grid_Y] = meshgrid(map.Grid_axis);
    map.Grid_r = sqrt(map.Grid_X.^2 + map.Grid_Y.^2);
    
    % recentering the map by a round number of pixel
    if  length(p.Results.shift) ~= 2
        error('Add_map(): for the centering option, a vector of 2 values must be given, for example[0.005 -0.02]')
    end
      
    map.offset_X = round(p.Results.shift(1)/map.res);
    map.offset_Y = round(p.Results.shift(2)/map.res);
    
    map.loaded = circshift(map.loaded,[-map.offset_Y map.offset_X]);
    % just shift by an integer number of pixel.     
    %figure(1); imagesc(map.Grid_axis,map.Grid_axis,map.loaded); axis square
    
    % If desired, remove the tilt
    if  ~isempty(p.Results.remove_tilt)
        diam = p.Results.remove_tilt;
        
        % Take all the points within the diameter
        map.central_ind =  intersect(find(map.Grid_r < diam/2),find(~isnan(map.loaded)));
        
        % Create the 1D vector for the fit
        map.fit_grid(:,1) = map.Grid_X(map.central_ind);
        map.fit_grid(:,2) = map.Grid_Y(map.central_ind);
        map.funcv = map.loaded(map.central_ind);
        
        % Function to fit
        func_curv = @(c,xdata)c(1)*xdata(:,1) + c(2)*xdata(:,2) + c(3);
        
        options = optimset('Display','off','MaxFunEvals',1E6,'TolFun',1E-12,'DiffMinChange',1E-12);
        
        c0 = [0 0 0];
        
        [map.fit_para,~,~,~,~] = lsqcurvefit(func_curv,c0,map.fit_grid,map.funcv,[],[],options);
        
        if p.Results.verbose
            fprintf('Substracted horizontal tilt [nrad]: %g \n',map.fit_para(1)*1E9)
            fprintf('Substracted vertical tilt [nrad]: %g \n',map.fit_para(2)*1E9)
        end
        
        map.fit = func_curv(map.fit_para,[map.Grid_X(:) map.Grid_Y(:)]);
        map.fit = reshape(map.fit,size(map.loaded));
        
        map.loaded = map.loaded - map.fit;
        
        %figure(2); imagesc(map.Grid_axis,map.Grid_axis,map.loaded); axis square
        
    end
    
    % If desired, remove the tilt and focus
    if  ~isempty(p.Results.remove_tilt_focus)
        diam = p.Results.remove_tilt_focus;
        
        % Take all the points within the diameter
        map.central_ind =  intersect(find(map.Grid_r < diam/2),find(~isnan(map.loaded)));
        
        % Create the 1D vector for the fit
        map.fit_grid(:,1) = map.Grid_X(map.central_ind);
        map.fit_grid(:,2) = map.Grid_Y(map.central_ind);
        map.funcv = map.loaded(map.central_ind);
        
        % Function to fit
        func_curv = @(c,xdata)c(1)*xdata(:,1) + c(2)*xdata(:,2) + c(3) + c(4)*(xdata(:,1).^2+xdata(:,2).^2);
        
        options = optimset('Display','off','MaxFunEvals',1E6,'TolFun',1E-12,'DiffMinChange',1E-12);
        
        c0 = [0 0 0 1/2000];
        
        [map.fit_para,~,~,~,~] = lsqcurvefit(func_curv,c0,map.fit_grid,map.funcv,[],[],options);
        
        if p.Results.verbose
            fprintf('Substracted radius of curvature [m]: %g \n',1/(2*map.fit_para(4)))
            fprintf('Substracted horizontal tilt [nrad]: %g \n',map.fit_para(1)*1E9)
            fprintf('Substracted vertical tilt [nrad]: %g \n',map.fit_para(2)*1E9)
        end
        
        map.fit = func_curv(map.fit_para,[map.Grid_X(:) map.Grid_Y(:)]);
        map.fit = reshape(map.fit,size(map.loaded));
        
        map.loaded = map.loaded - map.fit;
        
        %figure(2); imagesc(map.Grid_axis,map.Grid_axis,map.loaded); axis square
        
    end
    
    %figure(2); imagesc(map.Grid_axis,map.Grid_axis,map.loaded); axis square
    edge_value = find_edge_value(map.loaded);
    %edge_value = -1.81E-8;
    % Resample the loaded map to the grid of the interface
    map.resampled = interp2(map.Grid_X+p.Results.shift(1),map.Grid_Y-p.Results.shift(2),map.loaded,Iin.Grid.D2_X,Iin.Grid.D2_Y,'linear',edge_value);
    
elseif (n==2)    % The matrix is a 2 vector column, first column radius, second column sagitta change
    map.resampled = interp1(map.loaded(:,1),map.loaded(:,2),sqrt(Iin.Grid.D2_X.^2 + Iin.Grid.D2_Y.^2),'linear',0);
    map.resampled = map.resampled - map.loaded(1,2);
    
    % Rescale
    if  ~isempty(p.Results.scale)
        map.resampled = map.resampled * p.Results.scale;
    end
    edge_value = map.resampled(end);
else
    error('Add_map(): the loaded map is not square or a two columns array')
end

% remove the offset in the map:
%map.resampled = map.resampled - mean(map.resampled(~isnan(map.resampled)));

% Fix the RMS
if  ~isempty(p.Results.RMS)
    map.resampled = map.resampled/std(map.resampled(map.central_ind)) * p.Results.RMS;
end

map.resampled = rot90(map.resampled,round(p.Results.rotate));

% Remove the possible NaN remaining
map.resampled(isnan(map.resampled)) = edge_value;

Iout.surface =  Iin.surface - map.resampled;

end
