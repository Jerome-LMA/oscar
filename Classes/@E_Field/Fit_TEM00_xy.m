function varargout = Fit_TEM00_xy(Ein)
%  Fit_TEM00_xy() find the beam radius and the wavefront curvature for a TEM00
%  beam with astigmatism
%
% Required the 'Optimisation Toolbox' to work

if (~exist('lsqcurvefit','file'))
    error('Fit_TEM00() required the Optimisation Toolbox to run')
end

p  = inputParser;
p.FunctionName = 'Fit the TEM00 beam parameters';

% Check if the first argument is an E_Field
p.addRequired('Ein', @(x)isa(x, 'E_Field'));
p.parse(Ein)

Ein = p.Results.Ein;

if Calculate_power(Ein)==0
    error('Fit_TEM00_xy(): No power in the field')
end

% Find if the mode is a TEM00 or not
[~,m,n] = Read_mode_name(Ein.Mode_name);

if (m~=0) || (n~=0)
    error('Fit_TEM00_xy(): Can only be used to fit the fundamental mode. For higher order mode use Fit_E_Field()')
end

% For better accuracy, normalize the beam
Ein =  Normalise_E(Ein);

Fit_radius = Ein.Grid.Length/2;

% Take only the central part of the grid
tmp_index = find(Ein.Grid.D2_square < Fit_radius^2);
tmp_grid(:,1) = (Ein.Grid.D2_X(tmp_index)).^2;
tmp_grid(:,2) = (Ein.Grid.D2_Y(tmp_index)).^2;
Power_distri = abs(Ein.Field(tmp_index)).^2;

% Fit the power first
func_gauss = @(c,xdata) c(1)*(exp(-2*(xdata(:,1)/(c(2)^2) + xdata(:,2)/(c(3)^2)) ));

options = optimset('Display','off','MaxFunEvals',1E6,'TolFun',1E-8,'DiffMinChange',1E-8);
c0 = [max(max(Power_distri)) Fit_radius/2 Fit_radius/2];

[Map.fit_para,~,~,~,~] = lsqcurvefit(func_gauss,c0,tmp_grid,Power_distri,[],[],options);

Beam_rad_x = abs(Map.fit_para(2));
Beam_rad_y = abs(Map.fit_para(3));

% Now find the RofC of the complex wavefront
clear tmp_index tmp_grid
tmp_index = find(Ein.Grid.D2_square < (mean([Beam_rad_x Beam_rad_y]*3))^2);
tmp_grid(:,1) = (Ein.Grid.D2_X(tmp_index)).^2;
tmp_grid(:,2) = (Ein.Grid.D2_Y(tmp_index)).^2;
tmp_amp = Ein.Field(tmp_index);

% first find a starting wavefront RofC
% So take a cross section of the phase in one direction and then in the other
Cross_sec_phase_x = unwrap(angle(Ein.Field(Ein.Grid.Half_num_point,:)));
Cross_sec_phase_y = unwrap(angle(Ein.Field(:,Ein.Grid.Half_num_point)));

% Will do the fit on the beam diameter, find the index
Cross_sec_index_x = intersect(find (Ein.Grid.Axis < Beam_rad_x),find (Ein.Grid.Axis > -Beam_rad_x));
Cross_sec_index_y = intersect(find (Ein.Grid.Axis < Beam_rad_y),find (Ein.Grid.Axis > -Beam_rad_y));

poly = polyfit(Ein.Grid.Axis(Cross_sec_index_x),Cross_sec_phase_x(Cross_sec_index_x),2);
beam_radius_fit_x = -Ein.k_prop/(2*poly(1)); % first approximation in 1D


Cross_sec_phase_y = Cross_sec_phase_y';
poly = polyfit(Ein.Grid.Axis(Cross_sec_index_y),Cross_sec_phase_y(Cross_sec_index_y),2);
beam_radius_fit_y = -Ein.k_prop/(2*poly(1)); % first approximation in 1D

% Then fit the phase in 2D
func_gauss = @(c,xdata) c(1)*exp(-(xdata(:,1)/Beam_rad_x^2 + xdata(:,2)/Beam_rad_y^2)).*...
    exp(-1i*Ein.k_prop*(xdata(:,1)/(2*c(2)) + xdata(:,2)/(2*c(3)) ));

options = optimset('Display','off','MaxFunEvals',1E6,'TolFun',1E-25,'DiffMinChange',1E-25);
c0 = [max(max(Ein.Field)) beam_radius_fit_x beam_radius_fit_y];

[Map.fit_para,~,~,~,~] = lsqcurvefit(func_gauss,c0,tmp_grid,tmp_amp,[],[],options);
Beam_RofC_x = real(Map.fit_para(2));
Beam_RofC_y = real(Map.fit_para(3));


switch nargout
    case 0
        fprintf('Beam radius x [m]: %g  \t \t Wavefront curvature x [m]: %g  \n',Beam_rad_x,Beam_RofC_x)
        fprintf('Beam radius y [m]: %g  \t \t Wavefront curvature y [m]: %g  \n',Beam_rad_y,Beam_RofC_y)
    case 1
        Beam_rad(1,:) = Beam_rad_x;
        Beam_rad(2,:) = Beam_rad_y;
        varargout{1} = Beam_rad;
    case 2
        Beam_rad(1,:) = Beam_rad_x;
        Beam_rad(2,:) = Beam_rad_y;
        
        Beam_RofC(1,:) = Beam_RofC_x;
        Beam_RofC(2,:) = Beam_RofC_y;
        
        varargout{1} = Beam_rad;
        varargout{2} = Beam_RofC;
    otherwise
        error('Fit_TEM00(): Too many output argument')
end



end
